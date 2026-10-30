from __future__ import annotations

import csv
import hashlib
import io
import json
import re
import sqlite3
from contextlib import asynccontextmanager
from pathlib import Path
from typing import Any, Iterable

from fastapi import Body, FastAPI, HTTPException, Query, Request, Response
from fastapi.middleware.cors import CORSMiddleware
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel

from app.chemistry import Residue, aligned_pair_svgs, parent_scaffold_svg, residue_building_block_svg, structure_svg
from app.config import MMP_ENGINE_VERSION, PARSER_VERSION, STANDARDIZATION_VERSION, settings
from app.db import connect, initialize, row_dict, transaction, utcnow
from app.project_names import unique_project_name
from app.schemas import ImportResponse, ProjectCreate, ProjectUpdate
from app.sequence import building_block_token, residue_token, sequence_from_decomposition, sequence_with_site
from app.service import stage_import


@asynccontextmanager
async def lifespan(_app: FastAPI):
    initialize()
    yield


app = FastAPI(title="PepSeMan", version="0.1.0", lifespan=lifespan)
app.add_middleware(
    CORSMiddleware,
    allow_origins=["http://localhost:5173", "http://127.0.0.1:5173"],
    allow_credentials=False,
    allow_methods=["*"],
    allow_headers=["*"],
)


def require_project(db: sqlite3.Connection, project_id: int) -> sqlite3.Row:
    row = db.execute("SELECT * FROM projects WHERE project_id=?", (project_id,)).fetchone()
    if row is None:
        raise HTTPException(404, "Project not found")
    return row


def analysis_versions(db: sqlite3.Connection, project_id: int) -> tuple[str, str]:
    active = db.execute(
        """SELECT parser_version, mmp_engine_version FROM analysis_revisions
           WHERE project_id=? ORDER BY revision_id DESC LIMIT 1""",
        (project_id,),
    ).fetchone()
    return (active["parser_version"], active["mmp_engine_version"]) if active else (PARSER_VERSION, MMP_ENGINE_VERSION)


def sequence_aliases(db: sqlite3.Connection) -> dict[str, str]:
    return {row["variant_hash"]: row["abbreviation"] for row in db.execute("SELECT variant_hash, abbreviation FROM sequence_aliases")}


def _dictionary_rows(db: sqlite3.Connection) -> list[dict[str, Any]]:
    rows = db.execute(
        """WITH representatives AS (
               SELECT r.*, ROW_NUMBER() OVER (
                   PARTITION BY r.variant_hash ORDER BY r.parser_version DESC, r.structure_hash, r.position
               ) AS representative_rank
               FROM residues r WHERE r.parser_version IN (?, 'backbone-v1')
           )
           SELECT r.variant_hash, r.residue_type, r.r1, r.r2, r.r3, r.bridge,
                  r.stereo_status, a.abbreviation AS custom_abbreviation
           FROM representatives r
           LEFT JOIN sequence_aliases a ON a.variant_hash=r.variant_hash
           WHERE r.representative_rank=1 ORDER BY r.variant_hash""",
        (PARSER_VERSION,),
    ).fetchall()
    result = []
    for row in rows:
        item = dict(row)
        residue = Residue(0, item["residue_type"], item["r1"], item["r2"], item["r3"], item["bridge"], item["stereo_status"], [], {}, item["variant_hash"])
        item["default_abbreviation"] = residue_token(residue)
        if not (item["default_abbreviation"].startswith("[") and item["default_abbreviation"].endswith("]")):
            continue
        item["abbreviation"] = residue_token(residue, {item["variant_hash"]: item["custom_abbreviation"]}) if item["custom_abbreviation"] else item["default_abbreviation"]
        item["structure_url"] = f"/api/sequence-dictionary/{item['variant_hash']}/svg"
        result.append(item)
    return result


class SequenceAliasUpdate(BaseModel):
    abbreviation: str | None = None


@app.get("/api/sequence-dictionary")
def sequence_dictionary() -> list[dict[str, Any]]:
    with connect() as db:
        return _dictionary_rows(db)


@app.put("/api/sequence-dictionary/{variant_hash}")
def update_sequence_alias(variant_hash: str, payload: SequenceAliasUpdate) -> dict[str, Any]:
    if not re.fullmatch(r"[0-9a-f]{64}", variant_hash):
        raise HTTPException(422, "Invalid residue variant hash")
    alias = (payload.abbreviation or "").strip()
    if alias.startswith("[") and alias.endswith("]"):
        alias = alias[1:-1]
    if alias and not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_-]{0,15}", alias):
        raise HTTPException(422, "Use 1–16 ASCII letters, digits, _ or -; no spaces")
    with connect() as db:
        entries = _dictionary_rows(db)
        current = next((item for item in entries if item["variant_hash"] == variant_hash), None)
        if current is None:
            raise HTTPException(404, "Residue building block not found")
        if alias:
            if len(alias) == 1 and alias != current["default_abbreviation"]:
                raise HTTPException(422, "One-letter codes are reserved for the matching natural amino acid")
            if alias.casefold() in {"ac", "h", "oh", "nh2"}:
                raise HTTPException(422, "Terminal-cap codes are reserved")
            for item in entries:
                if item["variant_hash"] != variant_hash and item["default_abbreviation"].strip("[]").casefold() == alias.casefold():
                    raise HTTPException(409, "Abbreviation conflicts with an automatic building-block code")
        with transaction(db):
            try:
                if alias:
                    db.execute(
                        """INSERT INTO sequence_aliases(variant_hash, abbreviation, updated_at)
                           VALUES (?, ?, ?) ON CONFLICT(variant_hash) DO UPDATE SET
                           abbreviation=excluded.abbreviation, updated_at=excluded.updated_at""",
                        (variant_hash, alias, utcnow()),
                    )
                else:
                    db.execute("DELETE FROM sequence_aliases WHERE variant_hash=?", (variant_hash,))
            except sqlite3.IntegrityError as exc:
                raise HTTPException(409, "Abbreviation is already used by another building block") from exc
        return next(item for item in _dictionary_rows(db) if item["variant_hash"] == variant_hash)


@app.get("/api/sequence-dictionary/{variant_hash}/svg")
def sequence_dictionary_svg(variant_hash: str) -> Response:
    if not re.fullmatch(r"[0-9a-f]{64}", variant_hash):
        raise HTTPException(422, "Invalid residue variant hash")
    with connect() as db:
        row = db.execute(
            """SELECT r.position, d.decomposition_json, m.canonical_isomeric_smiles
               FROM residues r JOIN molecules m ON m.structure_hash=r.structure_hash
               JOIN decompositions d ON d.structure_hash=r.structure_hash
                   AND d.parser_version=r.parser_version AND d.standardization_version=?
               WHERE r.variant_hash=? AND r.parser_version IN (?, 'backbone-v1')
               ORDER BY r.parser_version DESC, r.structure_hash, r.position LIMIT 1""",
            (STANDARDIZATION_VERSION, variant_hash, PARSER_VERSION),
        ).fetchone()
        if row is None:
            raise HTTPException(404, "Residue building block not found")
        decomposition = json.loads(row["decomposition_json"])
        residue = Residue(**decomposition["residues"][row["position"] - 1])
        attachments = bridge_attachment_atoms(decomposition, row["position"])
        svg = residue_building_block_svg(row["canonical_isomeric_smiles"], residue, width=240, height=150, attachment_atoms=attachments)
        return Response(svg, media_type="image/svg+xml", headers={"Cache-Control": "private, max-age=3600"})


def parse_sort(sort: str | None, allowed: dict[str, str], default: str) -> str:
    if not sort:
        return default
    terms: list[str] = []
    for item in sort.split(","):
        name, separator, direction = item.partition(":")
        if name not in allowed or (separator and direction not in {"asc", "desc"}):
            raise HTTPException(422, f"Invalid sort term: {item}")
        expression = allowed[name]
        order = direction.upper() if separator else "ASC"
        terms.append(f"({expression} IS NULL) ASC")
        terms.append(f"{expression} {order}")
    return ", ".join(terms)


def sort_display_rows(rows: list[dict[str, Any]], terms: list[tuple[str, bool]], stable_key: str = "pair_key") -> list[dict[str, Any]]:
    """Sort computed MMP columns after reference orientation and sequence rendering."""
    rows.sort(key=lambda row: str(row[stable_key]))
    for name, descending in reversed(terms):
        present = [row for row in rows if row[name] is not None]
        missing = [row for row in rows if row[name] is None]
        present.sort(
            key=lambda row: (0, row[name]) if isinstance(row[name], (int, float)) else (1, str(row[name]).casefold()),
            reverse=descending,
        )
        rows = present + missing
    return rows


def page_query(
    db: sqlite3.Connection,
    sql: str,
    count_sql: str,
    parameters: Iterable[Any],
    order_by: str,
    page: int,
    page_size: int,
) -> dict[str, Any]:
    values = tuple(parameters)
    total = db.execute(count_sql, values).fetchone()[0]
    rows = db.execute(
        f"{sql} ORDER BY {order_by} LIMIT ? OFFSET ?",
        (*values, page_size, (page - 1) * page_size),
    ).fetchall()
    return {"rows": [dict(row) for row in rows], "total": total, "page": page, "page_size": page_size}


@app.get("/api/health")
def health() -> dict[str, str]:
    instance_id = hashlib.sha256(str(settings.database_path).encode()).hexdigest()[:16]
    return {"status": "ok", "instance_id": instance_id}


@app.get("/api/projects")
def projects(sort: str | None = None) -> list[dict[str, Any]]:
    allowed = {name: name for name in ("project_id", "name", "settings_json", "created_at", "updated_at")}
    order = parse_sort(sort, allowed, "updated_at DESC, project_id DESC")
    with connect() as db:
        return [dict(row) for row in db.execute(f"SELECT * FROM projects ORDER BY {order}")]


@app.post("/api/projects", status_code=201)
def create_project(payload: ProjectCreate) -> dict[str, Any]:
    now = utcnow()
    with connect() as db:
        name = unique_project_name(db, payload.name, now)
        cursor = db.execute(
            "INSERT INTO projects(name, settings_json, created_at, updated_at) VALUES (?, '{}', ?, ?)",
            (name, now, now),
        )
        return dict(db.execute("SELECT * FROM projects WHERE project_id=?", (cursor.lastrowid,)).fetchone())


@app.patch("/api/projects/{project_id}")
def rename_project(project_id: int, payload: ProjectUpdate) -> dict[str, Any]:
    with connect() as db:
        project = require_project(db, project_id)
        name = unique_project_name(db, payload.name, project["created_at"], exclude_id=project_id)
        db.execute(
            "UPDATE projects SET name=?, updated_at=? WHERE project_id=?",
            (name, utcnow(), project_id),
        )
        return dict(db.execute("SELECT * FROM projects WHERE project_id=?", (project_id,)).fetchone())


@app.delete("/api/projects/{project_id}", status_code=204)
def delete_project(project_id: int) -> Response:
    with connect() as db:
        require_project(db, project_id)
        db.execute("DELETE FROM projects WHERE project_id=?", (project_id,))
    return Response(status_code=204)


@app.post("/api/projects/{project_id}/imports", response_model=ImportResponse, status_code=202)
async def import_file(
    project_id: int,
    request: Request,
    source_name: str = Query(min_length=1, max_length=255),
    source_format: str = Query(pattern="^(csv|smi|sdf)$"),
) -> dict[str, Any]:
    content = await request.body()
    if not content:
        raise HTTPException(400, "Empty upload")
    if len(content) > 200 * 1024 * 1024:
        raise HTTPException(413, "Prototype upload limit is 200 MiB")
    with connect() as db:
        try:
            return stage_import(db, project_id, source_name, source_format, content)
        except LookupError as exc:
            raise HTTPException(404, str(exc)) from exc


COMMON_PAGE = {
    "page": Query(1, ge=1),
    "page_size": Query(50, ge=1, le=500),
}


@app.get("/api/projects/{project_id}/imports")
def list_imports(
    project_id: int,
    page: int = Query(1, ge=1),
    page_size: int = Query(50, ge=1, le=500),
    sort: str | None = None,
) -> dict[str, Any]:
    columns = {
        name: name
        for name in (
            "import_id",
            "source_name",
            "source_format",
            "file_sha256",
            "imported_at",
            "import_status",
            "dataset_revision",
            "summary_json",
        )
    }
    with connect() as db:
        require_project(db, project_id)
        return page_query(
            db,
            "SELECT import_id, source_name, source_format, file_sha256, imported_at, import_status, dataset_revision, summary_json FROM imports WHERE project_id=?",
            "SELECT count(*) FROM imports WHERE project_id=?",
            (project_id,),
            parse_sort(sort, columns, "dataset_revision DESC"),
            page,
            page_size,
        )


@app.get("/api/projects/{project_id}/qc")
def qc_rows(
    project_id: int,
    page: int = Query(1, ge=1),
    page_size: int = Query(50, ge=1, le=500),
    sort: str | None = None,
    status: str | None = None,
    id_filter: str | None = Query(None, max_length=120),
) -> dict[str, Any]:
    columns = {
        name: name
        for name in (
            "record_id",
            "source_molecule_id",
            "status",
            "exclusion_reason",
            "structure_hash",
            "original_smiles",
        )
    }
    params: list[Any] = [project_id]
    predicate = "project_id=?"
    if status:
        predicate += " AND status=?"
        params.append(status)
    if id_filter and id_filter.strip():
        predicate += " AND instr(lower(coalesce(source_molecule_id, '')), lower(?)) > 0"
        params.append(id_filter.strip())
    with connect() as db:
        require_project(db, project_id)
        return page_query(
            db,
            f"SELECT record_id, source_molecule_id, original_smiles, structure_hash, status, exclusion_reason FROM source_records WHERE {predicate}",
            f"SELECT count(*) FROM source_records WHERE {predicate}",
            params,
            parse_sort(sort, columns, "record_id ASC"),
            page,
            page_size,
        )


@app.get("/api/projects/{project_id}/decompositions")
def decompositions(
    project_id: int,
    page: int = Query(1, ge=1),
    page_size: int = Query(50, ge=1, le=500),
    sort: str | None = None,
    all_rows: bool = False,
    value_metric: str = Query("pchembl", pattern="^(pchembl|molecular_weight|length)$"),
    sequence_filter: str = Query("", max_length=200),
    cyclization_type: str = Query("all", pattern="^(all|linear|head_to_tail|thioether|disulfide|olefin_metathesis|triazole)$"),
) -> dict[str, Any]:
    value_expressions = {
        "pchembl": """(SELECT a.pchembl FROM activities a
                      WHERE a.project_id=pm.project_id AND a.molecule_id=pm.molecule_id
                        AND a.active=1 AND a.relation='=' AND a.pchembl IS NOT NULL
                      ORDER BY a.pchembl DESC LIMIT 1)""",
        "molecular_weight": "m.molecular_weight",
        "length": "d.length",
    }
    columns = {
        "source_molecule_id": "pm.source_molecule_id",
        "n_pairs": "n_pairs",
        "value": "value",
        "structure_hash": "m.structure_hash",
        "length": "d.length",
        "n_cap": "d.n_cap",
        "c_cap": "d.c_cap",
        "status": "d.status",
        "backbone_kind": "d.backbone_kind",
        "topology_signature": "d.topology_signature",
        "molecular_weight": "m.molecular_weight",
        "canonical_isomeric_smiles": "m.canonical_isomeric_smiles",
    }
    base = """FROM project_molecules pm
              JOIN molecules m ON m.molecule_id=pm.molecule_id
              JOIN decompositions d ON d.structure_hash=m.structure_hash AND d.parser_version=?
              WHERE pm.project_id=? AND pm.active=1"""
    with connect() as db:
        require_project(db, project_id)
        parser_version, engine_version = analysis_versions(db, project_id)
        params = (parser_version, project_id)
        computed_columns = {"sequence", "cyclization_types", "link_count"}
        sort_terms: list[tuple[str, bool]] = []
        for term in sort.split(",") if sort else []:
            name, separator, direction = term.partition(":")
            if (name not in columns and name not in computed_columns) or (separator and direction not in {"asc", "desc"}):
                raise HTTPException(422, f"Invalid sort term: {term}")
            sort_terms.append((name, direction == "desc"))
        order = "n_pairs DESC, pm.source_molecule_id ASC" if any(name in computed_columns for name, _ in sort_terms) else parse_sort(sort, columns, "n_pairs DESC, pm.source_molecule_id ASC")
        rows = db.execute(
            f"""SELECT pm.source_molecule_id,
                       (SELECT count(*) FROM pairs p
                        WHERE p.project_id=pm.project_id AND p.engine_version=?
                          AND (p.first_structure_hash=m.structure_hash OR p.second_structure_hash=m.structure_hash)) AS n_pairs,
                       {value_expressions[value_metric]} AS value,
                       m.structure_hash, m.canonical_isomeric_smiles, m.molecular_weight,
                       d.length, d.n_cap, d.c_cap, d.status, d.backbone_kind,
                       d.topology_signature, d.cyclization_links_json, d.decomposition_json
                {base} ORDER BY {order}""",
            (engine_version, *params),
        ).fetchall()
        aliases = sequence_aliases(db)
        result_rows = []
        for row in rows:
            item = dict(row)
            decomposition = json.loads(item.pop("decomposition_json"))
            links = json.loads(item.pop("cyclization_links_json"))
            result_rows.append({
                "source_molecule_id": item.pop("source_molecule_id"),
                "n_pairs": item.pop("n_pairs"),
                "value": item.pop("value"),
                "sequence": sequence_from_decomposition(decomposition, aliases),
                "canonical_start": decomposition.get("canonical_start", 1),
                "cyclization_types": ", ".join(sorted({link["kind"] for link in links})) or "linear",
                "link_count": len(links),
                "cyclization_links": links,
                **item,
            })
        if sequence_filter.strip():
            needle = sequence_filter.strip()
            result_rows = [item for item in result_rows if needle in item["sequence"]]
        if cyclization_type != "all":
            result_rows = [item for item in result_rows if
                           (not item["cyclization_links"] if cyclization_type == "linear" else
                            any(link["kind"] == cyclization_type for link in item["cyclization_links"]))]
        if any(name in computed_columns for name, _ in sort_terms):
            result_rows = sort_display_rows(result_rows, sort_terms, "structure_hash")
        total = len(result_rows)
        if not all_rows:
            result_rows = result_rows[(page - 1) * page_size:page * page_size]
        return {"rows": result_rows, "total": total, "page": 1 if all_rows else page, "page_size": max(total, 1) if all_rows else page_size}


PAIR_SELECT = """
    SELECT p.pair_key, p.first_structure_hash, p.second_structure_hash,
           p.site_kind, p.site_position, p.second_site_position, p.context_hash,
           (SELECT group_concat(DISTINCT pm.source_molecule_id) FROM project_molecules pm
            WHERE pm.project_id=p.project_id AND pm.molecule_id=fm.molecule_id AND pm.active=1) AS before_reference,
           (SELECT group_concat(DISTINCT pm.source_molecule_id) FROM project_molecules pm
            WHERE pm.project_id=p.project_id AND pm.molecule_id=sm.molecule_id AND pm.active=1) AS after_reference,
           fm.canonical_isomeric_smiles AS before_smiles,
           sm.canonical_isomeric_smiles AS after_smiles,
           fr.residue_type AS before_residue_type, sr.residue_type AS after_residue_type,
           fr.r1 AS before_r1, sr.r1 AS after_r1,
           fr.r2 AS before_r2, sr.r2 AS after_r2,
           fr.r3 AS before_r3, sr.r3 AS after_r3,
           fr.bridge AS before_bridge, sr.bridge AS after_bridge,
           fm.molecular_weight AS before_mw, sm.molecular_weight AS after_mw,
           sm.molecular_weight - fm.molecular_weight AS delta_mw,
           fa.assay_id, fa.endpoint, fa.pchembl AS before_activity,
           sa.pchembl AS after_activity,
           CASE WHEN fa.relation='=' AND sa.relation='=' AND fa.assay_id=sa.assay_id
                THEN sa.pchembl-fa.pchembl ELSE NULL END AS delta_activity,
           fa.provenance AS source,
           '/api/projects/' || p.project_id || '/pairs/' || p.pair_key || '/preview' AS scaffold_preview
    FROM pairs p
    JOIN molecules fm ON fm.structure_hash=p.first_structure_hash
    JOIN molecules sm ON sm.structure_hash=p.second_structure_hash
    LEFT JOIN residues fr ON fr.structure_hash=p.first_structure_hash
         AND fr.parser_version=(SELECT ar.parser_version FROM analysis_revisions ar
                                WHERE ar.project_id=p.project_id ORDER BY ar.revision_id DESC LIMIT 1)
         AND fr.position=p.site_position AND p.site_kind='residue'
    LEFT JOIN residues sr ON sr.structure_hash=p.second_structure_hash
         AND sr.parser_version=(SELECT ar.parser_version FROM analysis_revisions ar
                                WHERE ar.project_id=p.project_id ORDER BY ar.revision_id DESC LIMIT 1)
         AND sr.position=COALESCE(p.second_site_position,p.site_position) AND p.site_kind='residue'
    LEFT JOIN activities fa ON fa.activity_id=(
        SELECT ax.activity_id FROM activities ax
        WHERE ax.project_id=p.project_id AND ax.molecule_id=fm.molecule_id
          AND ax.active=1
        ORDER BY (ax.relation='=') DESC, ax.pchembl DESC LIMIT 1
    )
    LEFT JOIN activities sa ON sa.activity_id=(
        SELECT ax.activity_id FROM activities ax
        WHERE ax.project_id=p.project_id AND ax.molecule_id=sm.molecule_id
          AND ax.active=1
          AND (fa.assay_id IS NULL OR ax.assay_id=fa.assay_id)
          AND (fa.endpoint IS NULL OR ax.endpoint=fa.endpoint)
        ORDER BY (ax.relation='=') DESC, ax.pchembl DESC LIMIT 1
    )
    WHERE p.project_id=? AND p.engine_version=?
"""


@app.get("/api/projects/{project_id}/pairs")
def pairs(
    project_id: int,
    page: int = Query(1, ge=1),
    page_size: int = Query(50, ge=1, le=500),
    sort: str | None = None,
    reference_hash: str | None = None,
) -> dict[str, Any]:
    columns = {
        name: name
        for name in (
            "pair_key",
            "first_structure_hash",
            "second_structure_hash",
            "site_kind",
            "site_position",
            "second_site_position",
            "context_hash",
            "before_reference",
            "after_reference",
            "before_smiles",
            "after_smiles",
            "before_residue_type",
            "after_residue_type",
            "before_r1",
            "after_r1",
            "before_r2",
            "after_r2",
            "before_r3",
            "after_r3",
            "before_bridge",
            "after_bridge",
            "before_mw",
            "after_mw",
            "delta_mw",
            "assay_id",
            "endpoint",
            "before_activity",
            "after_activity",
            "delta_activity",
            "source",
            "scaffold_preview",
        )
    }
    computed_columns = {"parent_sequence", "original", "convert", "reference_structure_hash", "backbone_kind", "topology_signature", "cyclization_links", "cyclization_types", "residue_count", "partner_site_position"}
    sort_terms: list[tuple[str, bool]] = []
    for term in sort.split(",") if sort else []:
        name, separator, direction = term.partition(":")
        if (name not in columns and name not in computed_columns) or (separator and direction not in {"asc", "desc"}):
            raise HTTPException(422, f"Invalid sort term: {term}")
        sort_terms.append((name, direction == "desc"))
    sort_after_render = any(
        name in computed_columns or (reference_hash and name.startswith(("before_", "after_", "delta_")))
        for name, _ in sort_terms
    )
    with connect() as db:
        require_project(db, project_id)
        parser_version, engine_version = analysis_versions(db, project_id)
        if reference_hash and db.execute(
            """SELECT 1 FROM contexts WHERE project_id=? AND structure_hash=?
               AND engine_version=? LIMIT 1""",
            (project_id, reference_hash, engine_version),
        ).fetchone() is None:
            raise HTTPException(404, "Selected structure is not an accepted peptide in this project")
        filter_sql = " AND (p.first_structure_hash=? OR p.second_structure_hash=?)" if reference_hash else ""
        parameters: tuple[Any, ...] = (project_id, engine_version)
        if reference_hash:
            parameters += (reference_hash, reference_hash)
        if sort_after_render:
            rows = [dict(row) for row in db.execute(PAIR_SELECT + filter_sql, parameters).fetchall()]
            result = {"rows": rows, "total": len(rows), "page": page, "page_size": page_size}
        else:
            result = page_query(
                db,
                PAIR_SELECT + filter_sql,
                "SELECT count(*) FROM pairs p WHERE p.project_id=? AND p.engine_version=?" + filter_sql,
                parameters,
                parse_sort(sort, columns, "site_position ASC, pair_key ASC"),
                page,
                page_size,
            )
        aliases = sequence_aliases(db)
        structure_hashes = {
            structure_hash
            for row in result["rows"]
            for structure_hash in (row["first_structure_hash"], row["second_structure_hash"])
        }
        decomposition_by_hash = {}
        if structure_hashes:
            placeholders = ",".join("?" for _ in structure_hashes)
            decomposition_by_hash = {
                row["structure_hash"]: json.loads(row["decomposition_json"])
                for row in db.execute(
                    f"SELECT structure_hash, decomposition_json FROM decompositions WHERE parser_version=? AND structure_hash IN ({placeholders})",
                    (parser_version, *structure_hashes),
                )
            }
        for row in result["rows"]:
            first_position = row["site_position"]
            second_position = row["second_site_position"] or first_position
            if reference_hash and row["second_structure_hash"] == reference_hash:
                for field in ("reference", "smiles", "residue_type", "r1", "r2", "r3", "bridge", "mw", "activity"):
                    row[f"before_{field}"], row[f"after_{field}"] = row[f"after_{field}"], row[f"before_{field}"]
                for field in ("delta_mw", "delta_activity"):
                    if row[field] is not None:
                        row[field] = -row[field]
                row["site_position"] = second_position
                row["partner_site_position"] = first_position
            else:
                row["partner_site_position"] = second_position
            row["reference_structure_hash"] = reference_hash or row["first_structure_hash"]
            parent_hash = row["reference_structure_hash"]
            partner_hash = row["second_structure_hash"] if parent_hash == row["first_structure_hash"] else row["first_structure_hash"]
            parent = decomposition_by_hash.get(parent_hash)
            partner = decomposition_by_hash.get(partner_hash)
            row["parent_sequence"], row["parent_sequence_highlight_start"], row["parent_sequence_highlight_end"] = (
                sequence_with_site(parent, row["site_kind"], row["site_position"], aliases) if parent else ("", 0, 0)
            )
            row["original"] = building_block_token(parent, row["site_kind"], row["site_position"], aliases) if parent else ""
            row["convert"] = building_block_token(partner, row["site_kind"], row["partner_site_position"], aliases) if partner else ""
            row["backbone_kind"] = parent.get("backbone_kind", "linear") if parent else "linear"
            row["canonical_start"] = parent.get("canonical_start", 1) if parent else 1
            row["residue_count"] = len(parent.get("residues", [])) if parent else 0
            row["cyclization_links"] = parent.get("cyclization_links", []) if parent else []
            row["cyclization_types"] = ", ".join(sorted({link["kind"] for link in row["cyclization_links"]})) or "linear"
            row["topology_signature"] = parent.get("topology_signature", "") if parent else ""
        if sort_after_render:
            result["rows"] = sort_display_rows(result["rows"], sort_terms)[(page - 1) * page_size:page * page_size]
        return result


def _residues(db: sqlite3.Connection, structure_hash: str, parser_version: str) -> list[Residue]:
    row = db.execute(
        "SELECT decomposition_json FROM decompositions WHERE structure_hash=? AND parser_version=?",
        (structure_hash, parser_version),
    ).fetchone()
    if row is None:
        return []
    return [Residue(**item) for item in json.loads(row["decomposition_json"])["residues"]]


def bridge_attachment_atoms(topology: dict[str, Any], position: int) -> list[int]:
    return [
        int(end["atom_index"])
        for link in topology.get("cyclization_links", []) if link["kind"] != "head_to_tail"
        for end in (link["first"], link["second"])
        if end["kind"] == "residue" and end["position"] == position
    ]


@app.get("/api/projects/{project_id}/pairs/selection")
def pair_selection(project_id: int, reference_hash: str) -> list[dict[str, str]]:
    with connect() as db:
        require_project(db, project_id)
        _, engine_version = analysis_versions(db, project_id)
        if db.execute(
            "SELECT 1 FROM contexts WHERE project_id=? AND structure_hash=? AND engine_version=? LIMIT 1",
            (project_id, reference_hash, engine_version),
        ).fetchone() is None:
            raise HTTPException(404, "Selected structure is not an accepted peptide in this project")
        return [dict(row) for row in db.execute(
            """SELECT pair_key, first_structure_hash, second_structure_hash FROM pairs
               WHERE project_id=? AND engine_version=?
                 AND (first_structure_hash=? OR second_structure_hash=?)
               ORDER BY pair_key""",
            (project_id, engine_version, reference_hash, reference_hash),
        )]


@app.get("/api/projects/{project_id}/pairs/{pair_key}")
def pair_detail(project_id: int, pair_key: str, reference_hash: str | None = None) -> dict[str, Any]:
    with connect() as db:
        require_project(db, project_id)
        parser_version, engine_version = analysis_versions(db, project_id)
        pair = db.execute(
            "SELECT * FROM pairs WHERE project_id=? AND pair_key=? AND engine_version=?",
            (project_id, pair_key, engine_version),
        ).fetchone()
        if pair is None:
            raise HTTPException(404, "Pair not found in this project")
        if reference_hash and reference_hash not in (pair["first_structure_hash"], pair["second_structure_hash"]):
            raise HTTPException(404, "Pair does not contain the selected structure")
        first_mol = db.execute(
            "SELECT * FROM molecules WHERE structure_hash=?", (pair["first_structure_hash"],)
        ).fetchone()
        second_mol = db.execute(
            "SELECT * FROM molecules WHERE structure_hash=?", (pair["second_structure_hash"],)
        ).fetchone()
        first_residues = _residues(db, pair["first_structure_hash"], parser_version)
        second_residues = _residues(db, pair["second_structure_hash"], parser_version)
        property_row = db.execute(
            PAIR_SELECT + " AND p.pair_key=?",
            (project_id, engine_version, pair_key),
        ).fetchone()
        first_caps = db.execute("SELECT n_cap, c_cap, decomposition_json FROM decompositions WHERE structure_hash=? AND parser_version=?", (pair["first_structure_hash"], parser_version)).fetchone()
        second_caps = db.execute("SELECT n_cap, c_cap, decomposition_json FROM decompositions WHERE structure_hash=? AND parser_version=?", (pair["second_structure_hash"], parser_version)).fetchone()
        first_topology = json.loads(first_caps["decomposition_json"]) if first_caps else {}
        second_topology = json.loads(second_caps["decomposition_json"]) if second_caps else {}
        first_site = pair["site_position"]
        second_site = pair["second_site_position"] or first_site
        if pair["site_kind"] == "residue" and first_residues and second_residues and not first_topology.get("cyclization_links"):
            before_svg, after_svg, warning = aligned_pair_svgs(
                first_mol["canonical_isomeric_smiles"],
                second_mol["canonical_isomeric_smiles"],
                first_residues,
                second_residues,
                first_site,
            )
        else:
            before_svg = structure_svg(first_mol["canonical_isomeric_smiles"])
            after_svg = structure_svg(second_mol["canonical_isomeric_smiles"])
            warning = "Cyclic or terminal-cap pair uses full-structure depiction."
        result = {
            **dict(pair),
            "before": {**dict(first_mol), "n_cap": first_caps["n_cap"], "c_cap": first_caps["c_cap"], "residues": [vars(item) for item in first_residues], "svg": before_svg,
                       "bb_svg": residue_building_block_svg(first_mol["canonical_isomeric_smiles"], first_residues[first_site - 1], attachment_atoms=bridge_attachment_atoms(first_topology, first_site)) if pair["site_kind"] == "residue" else None},
            "after": {**dict(second_mol), "n_cap": second_caps["n_cap"], "c_cap": second_caps["c_cap"], "residues": [vars(item) for item in second_residues], "svg": after_svg,
                      "bb_svg": residue_building_block_svg(second_mol["canonical_isomeric_smiles"], second_residues[second_site - 1], attachment_atoms=bridge_attachment_atoms(second_topology, second_site)) if pair["site_kind"] == "residue" else None},
            "alignment_warning": warning,
            "properties": {
                "before_reference": property_row["before_reference"],
                "after_reference": property_row["after_reference"],
                "before_chembl_value": property_row["before_activity"],
                "after_chembl_value": property_row["after_activity"],
                "delta_chembl_value": property_row["delta_activity"],
                "before_mw": property_row["before_mw"],
                "after_mw": property_row["after_mw"],
                "delta_mw": property_row["delta_mw"],
                "assay_id": property_row["assay_id"],
                "endpoint": property_row["endpoint"],
            },
        }
        if reference_hash == pair["second_structure_hash"]:
            result["before"], result["after"] = result["after"], result["before"]
            for field in ("reference", "chembl_value", "mw"):
                before_key, after_key = f"before_{field}", f"after_{field}"
                result["properties"][before_key], result["properties"][after_key] = result["properties"][after_key], result["properties"][before_key]
            for field in ("delta_chembl_value", "delta_mw"):
                if result["properties"][field] is not None:
                    result["properties"][field] = -result["properties"][field]
        parent_site = second_site if reference_hash == pair["second_structure_hash"] else first_site
        result["site_position"] = parent_site
        result["partner_site_position"] = first_site if reference_hash == pair["second_structure_hash"] else second_site
        parent_topology = second_topology if reference_hash == pair["second_structure_hash"] else first_topology
        result["backbone_kind"] = parent_topology.get("backbone_kind", "linear")
        result["cyclization_links"] = parent_topology.get("cyclization_links", [])
        result["topology_signature"] = parent_topology.get("topology_signature", "")
        result["canonical_start"] = parent_topology.get("canonical_start", 1)
        result["scaffold_svg"] = parent_scaffold_svg(
            result["before"]["canonical_isomeric_smiles"],
            first_residues if result["before"]["structure_hash"] == pair["first_structure_hash"] else second_residues,
            changed_position=parent_site,
            site_kind=pair["site_kind"],
            cyclization_links=result["cyclization_links"],
        )
        result["parent_sequence"], result["parent_sequence_highlight_start"], result["parent_sequence_highlight_end"] = (
            sequence_with_site(parent_topology, pair["site_kind"], parent_site, sequence_aliases(db))
            if parent_topology else ("", 0, 0)
        )
        return result


@app.get("/api/projects/{project_id}/pairs/{pair_key}/preview")
def pair_preview(project_id: int, pair_key: str) -> Response:
    with connect() as db:
        require_project(db, project_id)
        parser_version, engine_version = analysis_versions(db, project_id)
        pair = db.execute(
            "SELECT * FROM pairs WHERE project_id=? AND pair_key=? AND engine_version=?",
            (project_id, pair_key, engine_version),
        ).fetchone()
        if pair is None:
            raise HTTPException(404, "Pair not found in this project")
        molecule = db.execute(
            "SELECT canonical_isomeric_smiles FROM molecules WHERE structure_hash=?",
            (pair["first_structure_hash"],),
        ).fetchone()
        highlight: list[int] = []
        if pair["site_kind"] == "residue":
            residues = _residues(db, pair["first_structure_hash"], parser_version)
            if 0 < pair["site_position"] <= len(residues):
                highlight = residues[pair["site_position"] - 1].atom_indices
        svg = structure_svg(molecule["canonical_isomeric_smiles"], highlight, width=260, height=120)
        return Response(svg, media_type="image/svg+xml")


@app.get("/api/projects/{project_id}/jobs")
def jobs(
    project_id: int,
    page: int = Query(1, ge=1),
    page_size: int = Query(50, ge=1, le=500),
    sort: str | None = None,
) -> dict[str, Any]:
    columns = {
        name: name
        for name in (
            "job_id",
            "project_id",
            "import_id",
            "job_type",
            "status",
            "progress",
            "error",
            "metrics_json",
            "created_at",
            "updated_at",
        )
    }
    with connect() as db:
        require_project(db, project_id)
        return page_query(
            db,
            "SELECT * FROM jobs WHERE project_id=?",
            "SELECT count(*) FROM jobs WHERE project_id=?",
            (project_id,),
            parse_sort(sort, columns, "created_at DESC"),
            page,
            page_size,
        )


@app.get("/api/projects/{project_id}/revisions")
def revisions(
    project_id: int,
    page: int = Query(1, ge=1),
    page_size: int = Query(50, ge=1, le=500),
    sort: str | None = None,
) -> dict[str, Any]:
    columns = {
        name: name
        for name in (
            "revision_id",
            "project_id",
            "input_revision",
            "parser_version",
            "mmp_engine_version",
            "summary_json",
            "created_at",
        )
    }
    with connect() as db:
        require_project(db, project_id)
        order = parse_sort(sort, columns, "created_at DESC")
        return page_query(
            db,
            "SELECT * FROM analysis_revisions WHERE project_id=?",
            "SELECT count(*) FROM analysis_revisions WHERE project_id=?",
            (project_id,),
            order,
            page,
            page_size,
        )


@app.post("/api/projects/{project_id}/revisions/{revision_id}/restore", response_model=ImportResponse, status_code=202)
def restore_revision(project_id: int, revision_id: int) -> dict[str, Any]:
    with connect() as db, transaction(db):
        require_project(db, project_id)
        revision = db.execute(
            "SELECT * FROM analysis_revisions WHERE project_id=? AND revision_id=?",
            (project_id, revision_id),
        ).fetchone()
        if revision is None:
            raise HTTPException(404, "Revision not found in this project")
        source = db.execute(
            """SELECT * FROM imports WHERE project_id=? AND dataset_revision=?
               ORDER BY import_id DESC LIMIT 1""",
            (project_id, revision["input_revision"]),
        ).fetchone()
        if source is None or not Path(source["source_path"]).is_file():
            raise HTTPException(409, "Revision source file is no longer available")
        next_revision = db.execute(
            "SELECT COALESCE(MAX(dataset_revision), 0) + 1 FROM imports WHERE project_id=?",
            (project_id,),
        ).fetchone()[0]
        now = utcnow()
        imported = db.execute(
            """INSERT INTO imports(
                   project_id, source_name, source_format, file_sha256, source_path,
                   imported_at, import_status, dataset_revision, summary_json
               ) VALUES (?, ?, ?, ?, ?, ?, 'queued', ?, '{}')""",
            (
                project_id,
                f"Restore revision {revision_id} · {source['source_name']}",
                source["source_format"],
                source["file_sha256"],
                source["source_path"],
                now,
                next_revision,
            ),
        )
        job = db.execute(
            """INSERT INTO jobs(project_id, import_id, job_type, status, progress, created_at, updated_at)
               VALUES (?, ?, 'import_and_analyze', 'queued', 0, ?, ?)""",
            (project_id, imported.lastrowid, now, now),
        )
        return {
            "import_id": imported.lastrowid,
            "job_id": job.lastrowid,
            "duplicate": False,
            "dataset_revision": next_revision,
        }


class ExportRequest(BaseModel):
    structure_hashes: list[str] = []
    mode: str = "measurement_rows"


@app.post("/api/projects/{project_id}/export")
def export_csv(project_id: int, payload: ExportRequest) -> Response:
    if payload.mode not in {"measurement_rows", "unique_compounds", "selected_compounds"}:
        raise HTTPException(422, "Unsupported export mode")
    if payload.mode == "selected_compounds" and not payload.structure_hashes:
        raise HTTPException(422, "Selected compounds export requires structure hashes")
    with connect() as db:
        require_project(db, project_id)
        parser_version, _ = analysis_versions(db, project_id)
        aliases = sequence_aliases(db)
        params: list[Any] = [project_id]
        structure_filter = ""
        if payload.structure_hashes:
            placeholders = ",".join("?" for _ in payload.structure_hashes)
            structure_filter = f" AND m.structure_hash IN ({placeholders})"
            params.extend(payload.structure_hashes)
        if payload.mode == "selected_compounds":
            rows = db.execute(
                f"""SELECT pm.source_molecule_id, m.structure_hash, m.canonical_isomeric_smiles,
                           m.molecular_weight, d.length, d.n_cap, d.c_cap, d.status, d.decomposition_json,
                           d.backbone_kind, d.topology_signature, d.cyclization_links_json,
                           a.assay_id, a.endpoint, a.relation, a.value, a.unit, a.pchembl, a.provenance
                    FROM project_molecules pm JOIN molecules m ON m.molecule_id=pm.molecule_id
                    LEFT JOIN decompositions d ON d.structure_hash=m.structure_hash AND d.parser_version=?
                    LEFT JOIN activities a ON a.project_id=pm.project_id AND a.molecule_id=pm.molecule_id AND a.active=1
                    WHERE pm.project_id=? AND pm.active=1 {structure_filter}
                    ORDER BY pm.source_molecule_id, a.assay_id, a.activity_id""",
                (parser_version, *params),
            ).fetchall()
        elif payload.mode == "unique_compounds":
            rows = db.execute(
                f"""SELECT DISTINCT m.structure_hash, pm.source_molecule_id,
                           m.canonical_isomeric_smiles, m.molecular_weight,
                           d.length, d.n_cap, d.c_cap, d.status
                    FROM project_molecules pm JOIN molecules m ON m.molecule_id=pm.molecule_id
                    LEFT JOIN decompositions d ON d.structure_hash=m.structure_hash AND d.parser_version=?
                    WHERE pm.project_id=? AND pm.active=1 {structure_filter}
                    ORDER BY pm.source_molecule_id, m.structure_hash""",
                (parser_version, *params),
            ).fetchall()
        else:
            rows = db.execute(
                f"""SELECT m.structure_hash, sr.source_molecule_id, m.canonical_isomeric_smiles,
                           a.assay_id, a.endpoint, a.relation, a.relation_raw, a.value,
                           a.unit, a.pchembl, a.provenance
                    FROM activities a JOIN molecules m ON m.molecule_id=a.molecule_id
                    JOIN source_records sr ON sr.record_id=a.source_record_id
                    WHERE a.project_id=? AND a.active=1 {structure_filter}
                    ORDER BY sr.source_molecule_id, a.assay_id, a.activity_id""",
                params,
            ).fetchall()
    output = io.StringIO()
    if rows:
        export_rows = []
        for row in rows:
            item = dict(row)
            if payload.mode == "selected_compounds":
                decomposition_json = item.pop("decomposition_json")
                item["sequence"] = sequence_from_decomposition(json.loads(decomposition_json), aliases) if decomposition_json else ""
                item["cyclization_links_json"] = item["cyclization_links_json"] or "[]"
            export_rows.append(item)
        writer = csv.DictWriter(output, fieldnames=export_rows[0].keys())
        writer.writeheader()
        writer.writerows(export_rows)
    return Response(
        output.getvalue(),
        media_type="text/csv; charset=utf-8",
        headers={"Content-Disposition": f'attachment; filename="project-{project_id}-{payload.mode}.csv"'},
    )


frontend_dist = Path(__file__).parents[1] / "frontend" / "dist"
if frontend_dist.is_dir():
    app.mount("/", StaticFiles(directory=frontend_dist, html=True), name="frontend")

