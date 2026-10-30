from __future__ import annotations

import hashlib
import json
import shutil
import sqlite3
import time
from collections import Counter
from pathlib import Path
from typing import Any

from app.chemistry import Decomposition, decompose_smiles
from app.config import MMP_ENGINE_VERSION, PARSER_VERSION, STANDARDIZATION_VERSION, settings
from app.db import json_text, transaction, utcnow
from app.importers import SourceRecord, iter_records
from app.mmp import update_project_contexts


def _optional_float(value: str | None) -> float | None:
    if value in (None, "", "None"):
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def normalize_relation(value: str | None) -> str | None:
    if value is None:
        return None
    normalized = value.strip().strip('"').strip("'").strip()
    return normalized if normalized in {"=", "<", ">", "<=", ">="} else normalized or None


def _record_hash(record: SourceRecord) -> str:
    payload = json.dumps(record.properties, sort_keys=True, ensure_ascii=False, separators=(",", ":"))
    return hashlib.sha256(f"{record.row_number}\0{payload}".encode()).hexdigest()


def _load_cached(db: sqlite3.Connection, structure_hash: str) -> Decomposition | None:
    row = db.execute(
        """SELECT decomposition_json FROM decompositions
           WHERE structure_hash=? AND parser_version=? AND standardization_version=?""",
        (structure_hash, PARSER_VERSION, STANDARDIZATION_VERSION),
    ).fetchone()
    if row is None:
        return None
    data = json.loads(row["decomposition_json"])
    from app.chemistry import Residue

    data["residues"] = [Residue(**item) for item in data["residues"]]
    return Decomposition(**data)


def _store_decomposition(db: sqlite3.Connection, result: Decomposition) -> None:
    if result.structure_hash is None:
        return
    db.execute(
        """INSERT OR REPLACE INTO decompositions(
               structure_hash, parser_version, standardization_version, length,
               n_cap, c_cap, status, error, decomposition_json, updated_at,
               backbone_kind, topology_signature, cyclization_links_json
           ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)""",
        (
            result.structure_hash,
            PARSER_VERSION,
            STANDARDIZATION_VERSION,
            result.length,
            result.n_cap,
            result.c_cap,
            result.status,
            result.error,
            json_text(result.as_json()),
            utcnow(),
            result.backbone_kind,
            result.topology_signature,
            json_text(result.cyclization_links),
        ),
    )
    db.execute(
        "DELETE FROM cyclization_links WHERE structure_hash=? AND parser_version=?",
        (result.structure_hash, PARSER_VERSION),
    )
    db.executemany(
        """INSERT INTO cyclization_links(
               structure_hash, parser_version, link_index, kind,
               first_site_kind, first_position, second_site_kind, second_position,
               bond_order, stereo, linker_smiles, link_json
           ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)""",
        [
            (
                result.structure_hash, PARSER_VERSION, index, link["kind"],
                link["first"]["kind"], link["first"]["position"],
                link["second"]["kind"], link["second"]["position"],
                link["bond_order"], link["stereo"], link["linker_smiles"], json_text(link),
            )
            for index, link in enumerate(result.cyclization_links)
        ],
    )
    db.execute(
        "DELETE FROM residues WHERE structure_hash=? AND parser_version=?",
        (result.structure_hash, PARSER_VERSION),
    )
    db.executemany(
        """INSERT INTO residues(
               structure_hash, parser_version, position, residue_type, r1, r2,
               r3, bridge, stereo_status, atom_indices_json, variant_hash
           ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)""",
        [
            (
                result.structure_hash,
                PARSER_VERSION,
                residue.position,
                residue.residue_type,
                residue.r1,
                residue.r2,
                residue.r3,
                residue.bridge,
                residue.stereo_status,
                json_text(residue.atom_indices),
                residue.variant_hash,
            )
            for residue in result.residues
        ],
    )


def stage_import(
    db: sqlite3.Connection,
    project_id: int,
    source_name: str,
    source_format: str,
    content: bytes,
) -> dict[str, int | bool | None]:
    project = db.execute("SELECT 1 FROM projects WHERE project_id=?", (project_id,)).fetchone()
    if project is None:
        raise LookupError("Project not found")
    digest = hashlib.sha256(content).hexdigest()
    previous = db.execute(
        """SELECT * FROM imports WHERE project_id=? AND file_sha256=?
           ORDER BY dataset_revision DESC LIMIT 1""",
        (project_id, digest),
    ).fetchone()
    revision = db.execute(
        "SELECT COALESCE(MAX(dataset_revision), 0) + 1 AS value FROM imports WHERE project_id=?",
        (project_id,),
    ).fetchone()["value"]
    import_dir = settings.data_dir / "imports" / str(project_id)
    import_dir.mkdir(parents=True, exist_ok=True)
    safe_suffix = {"csv": ".csv", "smi": ".smi", "sdf": ".sdf"}[source_format]
    destination = import_dir / f"{digest}{safe_suffix}"
    if not destination.exists():
        destination.write_bytes(content)
    now = utcnow()
    with transaction(db):
        cursor = db.execute(
            """INSERT INTO imports(
                   project_id, source_name, source_format, file_sha256, source_path,
                   imported_at, import_status, dataset_revision, summary_json
               ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)""",
            (
                project_id,
                source_name,
                source_format,
                digest,
                str(destination),
                now,
                "reused" if previous else "queued",
                revision,
                previous["summary_json"] if previous else "{}",
            ),
        )
        import_id = cursor.lastrowid
        if previous:
            return {"import_id": import_id, "job_id": None, "duplicate": True, "dataset_revision": revision}
        job = db.execute(
            """INSERT INTO jobs(project_id, import_id, job_type, status, progress, created_at, updated_at)
               VALUES (?, ?, 'import_and_analyze', 'queued', 0, ?, ?)""",
            (project_id, import_id, now, now),
        )
    return {"import_id": import_id, "job_id": job.lastrowid, "duplicate": False, "dataset_revision": revision}


def process_import(db: sqlite3.Connection, import_id: int, progress=None) -> dict[str, Any]:
    started = time.perf_counter()
    import_row = db.execute("SELECT * FROM imports WHERE import_id=?", (import_id,)).fetchone()
    if import_row is None:
        raise LookupError("Import not found")
    project_id = import_row["project_id"]
    records = list(iter_records(Path(import_row["source_path"]), import_row["source_format"]))
    prepared: list[tuple[SourceRecord, Decomposition]] = []
    for index, record in enumerate(records, start=1):
        prepared.append((record, decompose_smiles(record.smiles)))
        if progress and (index % 20 == 0 or index == len(records)):
            # Chemistry is computed before the atomic database update so a
            # separate worker connection can publish durable progress without
            # contending with the import transaction's write lock.
            progress(0.05 + 0.55 * index / max(len(records), 1))
    status_counts: Counter[str] = Counter()
    accepted_hashes: set[str] = set()
    recomputed = 0
    reused = 0
    parseable = 0
    unique_source_ids: set[str] = set()

    with transaction(db):
        # A new non-duplicate import is a full dataset snapshot. Previous rows
        # remain available for provenance/revision audit, but current views use
        # only the atomically activated snapshot.
        db.execute("UPDATE activities SET active=0 WHERE project_id=?", (project_id,))
        db.execute("UPDATE project_molecules SET active=0 WHERE project_id=?", (project_id,))
        for index, (record, result) in enumerate(prepared, start=1):
            if result.structure_hash:
                parseable += 1
                cached = _load_cached(db, result.structure_hash)
                if cached is None:
                    _store_decomposition(db, result)
                    recomputed += 1
                else:
                    result = cached
                    reused += 1
            status_counts[result.status] += 1
            if record.source_molecule_id:
                unique_source_ids.add(record.source_molecule_id)
            cursor = db.execute(
                """INSERT INTO source_records(
                       project_id, import_id, original_row_data, original_row_hash,
                       source_molecule_id, source_activity_id, original_smiles,
                       structure_hash, status, exclusion_reason
                   ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)""",
                (
                    project_id,
                    import_id,
                    json_text(record.properties),
                    _record_hash(record),
                    record.source_molecule_id,
                    record.source_activity_id,
                    record.smiles,
                    result.structure_hash,
                    result.status,
                    result.error if result.status != "accepted" else None,
                ),
            )
            source_record_id = cursor.lastrowid
            if result.structure_hash and result.canonical_smiles:
                db.execute(
                    """INSERT OR IGNORE INTO molecules(
                           canonical_isomeric_smiles, structure_hash, source_metadata, molecular_weight
                       ) VALUES (?, ?, ?, ?)""",
                    (result.canonical_smiles, result.structure_hash, "{}", result.molecular_weight),
                )
                molecule_id = db.execute(
                    "SELECT molecule_id FROM molecules WHERE structure_hash=?", (result.structure_hash,)
                ).fetchone()["molecule_id"]
                db.execute(
                    """INSERT OR REPLACE INTO project_molecules(
                           project_id, molecule_id, source_molecule_id, active
                       ) VALUES (?, ?, ?, 1)""",
                    (project_id, molecule_id, record.source_molecule_id or ""),
                )
                if result.status == "accepted":
                    accepted_hashes.add(result.structure_hash)
                properties = record.properties
                db.execute(
                    """INSERT INTO activities(
                           project_id, molecule_id, assay_id, endpoint, relation, relation_raw,
                           value, unit, pchembl, provenance, source_record_id, active
                       ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, 1)""",
                    (
                        project_id,
                        molecule_id,
                        properties.get("Assay ChEMBL ID"),
                        properties.get("Standard Type"),
                        normalize_relation(properties.get("Standard Relation")),
                        properties.get("Standard Relation"),
                        _optional_float(properties.get("Standard Value")),
                        properties.get("Standard Units"),
                        _optional_float(properties.get("pChEMBL Value")),
                        json_text(
                            {
                                "document": properties.get("Document ChEMBL ID"),
                                "source": properties.get("Source Name"),
                                "source_name": import_row["source_name"],
                            }
                        ),
                        source_record_id,
                    ),
                )
        current_molecule_ids = {
            row["molecule_id"]
            for row in db.execute(
                """SELECT DISTINCT m.molecule_id FROM molecules m
                   JOIN source_records sr ON sr.structure_hash=m.structure_hash
                   WHERE sr.import_id=?""",
                (import_id,),
            )
        }
        context_metrics = update_project_contexts(db, project_id, accepted_hashes)
        elapsed = time.perf_counter() - started
        summary = {
            "records": len(records),
            "parseable": parseable,
            "source_molecule_ids": len(unique_source_ids),
            "unique_structures": len(current_molecule_ids),
            "accepted_structures": len(accepted_hashes),
            "status_counts": dict(sorted(status_counts.items())),
            "recomputed": recomputed,
            "reused": reused,
            "cache_reuse_rate": reused / max(reused + recomputed, 1),
            "elapsed_seconds": round(elapsed, 3),
            **context_metrics,
        }
        db.execute(
            "UPDATE imports SET import_status='complete', summary_json=? WHERE import_id=?",
            (json_text(summary), import_id),
        )
        db.execute(
            """INSERT INTO analysis_revisions(
                   project_id, input_revision, parser_version, mmp_engine_version,
                   summary_json, created_at
               ) VALUES (?, ?, ?, ?, ?, ?)""",
            (
                project_id,
                import_row["dataset_revision"],
                PARSER_VERSION,
                MMP_ENGINE_VERSION,
                json_text(summary),
                utcnow(),
            ),
        )
    return summary

