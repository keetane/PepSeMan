from __future__ import annotations

import json
import sqlite3
from contextlib import contextmanager
from datetime import UTC, datetime
from pathlib import Path
from typing import Any, Iterator

from app.config import MMP_ENGINE_VERSION, PARSER_VERSION, settings
from app.project_names import unique_project_name


SCHEMA_VERSION = 6


def utcnow() -> str:
    return datetime.now(UTC).isoformat()


def connect(path: Path | None = None) -> sqlite3.Connection:
    target = path or settings.database_path
    target.parent.mkdir(parents=True, exist_ok=True)
    connection = sqlite3.connect(target, timeout=30, isolation_level=None)
    connection.row_factory = sqlite3.Row
    connection.execute("PRAGMA journal_mode=WAL")
    connection.execute("PRAGMA foreign_keys=ON")
    connection.execute("PRAGMA busy_timeout=30000")
    return connection


@contextmanager
def transaction(connection: sqlite3.Connection) -> Iterator[sqlite3.Connection]:
    connection.execute("BEGIN IMMEDIATE")
    try:
        yield connection
    except Exception:
        connection.execute("ROLLBACK")
        raise
    else:
        connection.execute("COMMIT")


def json_text(value: Any) -> str:
    return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":"))


def initialize(path: Path | None = None) -> None:
    with connect(path) as db:
        db.executescript(
                """
                CREATE TABLE IF NOT EXISTS schema_migrations (
                    version INTEGER PRIMARY KEY,
                    applied_at TEXT NOT NULL
                );

                CREATE TABLE IF NOT EXISTS projects (
                    project_id INTEGER PRIMARY KEY AUTOINCREMENT,
                    name TEXT NOT NULL CHECK(length(trim(name)) > 0),
                    settings_json TEXT NOT NULL DEFAULT '{}',
                    created_at TEXT NOT NULL,
                    updated_at TEXT NOT NULL
                );

                CREATE TABLE IF NOT EXISTS imports (
                    import_id INTEGER PRIMARY KEY AUTOINCREMENT,
                    project_id INTEGER NOT NULL REFERENCES projects(project_id) ON DELETE CASCADE,
                    source_name TEXT NOT NULL,
                    source_format TEXT NOT NULL CHECK(source_format IN ('csv','smi','sdf')),
                    file_sha256 TEXT NOT NULL,
                    source_path TEXT NOT NULL,
                    imported_at TEXT NOT NULL,
                    import_status TEXT NOT NULL,
                    dataset_revision INTEGER NOT NULL,
                    summary_json TEXT NOT NULL DEFAULT '{}',
                    UNIQUE(project_id, file_sha256, dataset_revision)
                );

                CREATE TABLE IF NOT EXISTS source_records (
                    record_id INTEGER PRIMARY KEY AUTOINCREMENT,
                    project_id INTEGER NOT NULL REFERENCES projects(project_id) ON DELETE CASCADE,
                    import_id INTEGER NOT NULL REFERENCES imports(import_id) ON DELETE CASCADE,
                    original_row_data TEXT NOT NULL,
                    original_row_hash TEXT NOT NULL,
                    source_molecule_id TEXT,
                    source_activity_id TEXT,
                    original_smiles TEXT,
                    structure_hash TEXT,
                    status TEXT NOT NULL,
                    exclusion_reason TEXT,
                    UNIQUE(import_id, original_row_hash)
                );

                CREATE TABLE IF NOT EXISTS molecules (
                    molecule_id INTEGER PRIMARY KEY AUTOINCREMENT,
                    canonical_isomeric_smiles TEXT NOT NULL,
                    structure_hash TEXT NOT NULL UNIQUE,
                    source_metadata TEXT NOT NULL DEFAULT '{}',
                    molecular_weight REAL
                );

                CREATE TABLE IF NOT EXISTS project_molecules (
                    project_id INTEGER NOT NULL REFERENCES projects(project_id) ON DELETE CASCADE,
                    molecule_id INTEGER NOT NULL REFERENCES molecules(molecule_id) ON DELETE CASCADE,
                    source_molecule_id TEXT,
                    active INTEGER NOT NULL DEFAULT 1,
                    PRIMARY KEY(project_id, molecule_id, source_molecule_id)
                );

                CREATE TABLE IF NOT EXISTS activities (
                    activity_id INTEGER PRIMARY KEY AUTOINCREMENT,
                    project_id INTEGER NOT NULL REFERENCES projects(project_id) ON DELETE CASCADE,
                    molecule_id INTEGER NOT NULL REFERENCES molecules(molecule_id) ON DELETE CASCADE,
                    assay_id TEXT,
                    endpoint TEXT,
                    relation TEXT,
                    relation_raw TEXT,
                    value REAL,
                    unit TEXT,
                    pchembl REAL,
                    provenance TEXT NOT NULL,
                    source_record_id INTEGER NOT NULL REFERENCES source_records(record_id) ON DELETE CASCADE,
                    active INTEGER NOT NULL DEFAULT 1,
                    UNIQUE(project_id, source_record_id)
                );

                CREATE TABLE IF NOT EXISTS decompositions (
                    structure_hash TEXT NOT NULL,
                    parser_version TEXT NOT NULL,
                    standardization_version TEXT NOT NULL,
                    length INTEGER,
                    n_cap TEXT,
                    c_cap TEXT,
                    status TEXT NOT NULL,
                    error TEXT,
                    decomposition_json TEXT NOT NULL,
                    updated_at TEXT NOT NULL,
                    PRIMARY KEY(structure_hash, parser_version, standardization_version)
                );

                CREATE TABLE IF NOT EXISTS residues (
                    structure_hash TEXT NOT NULL,
                    parser_version TEXT NOT NULL,
                    position INTEGER NOT NULL,
                    residue_type TEXT NOT NULL,
                    r1 TEXT,
                    r2 TEXT,
                    r3 TEXT,
                    bridge TEXT,
                    stereo_status TEXT NOT NULL,
                    atom_indices_json TEXT NOT NULL,
                    variant_hash TEXT NOT NULL,
                    PRIMARY KEY(structure_hash, parser_version, position)
                );

                CREATE TABLE IF NOT EXISTS sequence_aliases (
                    variant_hash TEXT PRIMARY KEY,
                    abbreviation TEXT NOT NULL UNIQUE COLLATE NOCASE,
                    updated_at TEXT NOT NULL
                );
                CREATE INDEX IF NOT EXISTS idx_residue_variant ON residues(variant_hash, parser_version);

                CREATE TABLE IF NOT EXISTS contexts (
                    project_id INTEGER NOT NULL REFERENCES projects(project_id) ON DELETE CASCADE,
                    structure_hash TEXT NOT NULL,
                    site_kind TEXT NOT NULL,
                    site_position INTEGER NOT NULL,
                    context_hash TEXT NOT NULL,
                    variant_hash TEXT NOT NULL,
                    engine_version TEXT NOT NULL,
                    PRIMARY KEY(project_id, structure_hash, site_kind, site_position, engine_version)
                );
                CREATE INDEX IF NOT EXISTS idx_context_bucket
                    ON contexts(project_id, context_hash, engine_version);

                CREATE TABLE IF NOT EXISTS pairs (
                    project_id INTEGER NOT NULL REFERENCES projects(project_id) ON DELETE CASCADE,
                    pair_key TEXT NOT NULL,
                    first_structure_hash TEXT NOT NULL,
                    second_structure_hash TEXT NOT NULL,
                    site_kind TEXT NOT NULL,
                    site_position INTEGER NOT NULL,
                    context_hash TEXT NOT NULL,
                    engine_version TEXT NOT NULL,
                    created_at TEXT NOT NULL,
                    PRIMARY KEY(project_id, pair_key, engine_version),
                    CHECK(first_structure_hash < second_structure_hash)
                );

                CREATE TABLE IF NOT EXISTS jobs (
                    job_id INTEGER PRIMARY KEY AUTOINCREMENT,
                    project_id INTEGER NOT NULL REFERENCES projects(project_id) ON DELETE CASCADE,
                    import_id INTEGER REFERENCES imports(import_id) ON DELETE CASCADE,
                    job_type TEXT NOT NULL,
                    status TEXT NOT NULL,
                    progress REAL NOT NULL DEFAULT 0 CHECK(progress >= 0 AND progress <= 1),
                    error TEXT,
                    metrics_json TEXT NOT NULL DEFAULT '{}',
                    created_at TEXT NOT NULL,
                    updated_at TEXT NOT NULL
                );
                CREATE INDEX IF NOT EXISTS idx_jobs_claim ON jobs(status, created_at);

                CREATE TABLE IF NOT EXISTS analysis_revisions (
                    revision_id INTEGER PRIMARY KEY AUTOINCREMENT,
                    project_id INTEGER NOT NULL REFERENCES projects(project_id) ON DELETE CASCADE,
                    input_revision INTEGER NOT NULL,
                    parser_version TEXT NOT NULL,
                    mmp_engine_version TEXT NOT NULL,
                    summary_json TEXT NOT NULL,
                    created_at TEXT NOT NULL
                );

                CREATE TABLE IF NOT EXISTS depiction_cache (
                    pair_scaffold_hash TEXT NOT NULL,
                    renderer_version TEXT NOT NULL,
                    asset_path TEXT NOT NULL,
                    created_at TEXT NOT NULL,
                    PRIMARY KEY(pair_scaffold_hash, renderer_version)
                );

                CREATE TABLE IF NOT EXISTS cyclization_links (
                    structure_hash TEXT NOT NULL,
                    parser_version TEXT NOT NULL,
                    link_index INTEGER NOT NULL,
                    kind TEXT NOT NULL,
                    first_site_kind TEXT NOT NULL,
                    first_position INTEGER NOT NULL,
                    second_site_kind TEXT NOT NULL,
                    second_position INTEGER NOT NULL,
                    bond_order REAL NOT NULL,
                    stereo TEXT NOT NULL,
                    linker_smiles TEXT NOT NULL,
                    link_json TEXT NOT NULL,
                    PRIMARY KEY(structure_hash, parser_version, link_index)
                );
                CREATE INDEX IF NOT EXISTS idx_cyclization_kind ON cyclization_links(kind, parser_version);
                """
        )
        decomposition_columns = {row["name"] for row in db.execute("PRAGMA table_info(decompositions)")}
        for name, definition in (
            ("backbone_kind", "TEXT NOT NULL DEFAULT 'linear'"),
            ("topology_signature", "TEXT NOT NULL DEFAULT ''"),
            ("cyclization_links_json", "TEXT NOT NULL DEFAULT '[]'"),
        ):
            if name not in decomposition_columns:
                db.execute(f"ALTER TABLE decompositions ADD COLUMN {name} {definition}")
        pair_columns = {row["name"] for row in db.execute("PRAGMA table_info(pairs)")}
        if "second_site_position" not in pair_columns:
            db.execute("ALTER TABLE pairs ADD COLUMN second_site_position INTEGER")
        activity_columns = {row["name"] for row in db.execute("PRAGMA table_info(activities)")}
        if "active" not in activity_columns:
            db.execute("ALTER TABLE activities ADD COLUMN active INTEGER NOT NULL DEFAULT 1")
        if not db.execute("SELECT 1 FROM schema_migrations WHERE version=4").fetchone():
            with transaction(db):
                for project in db.execute("SELECT project_id, name, created_at FROM projects ORDER BY project_id").fetchall():
                    source = db.execute(
                        "SELECT source_name FROM imports WHERE project_id=? ORDER BY imported_at, import_id LIMIT 1",
                        (project["project_id"],),
                    ).fetchone()
                    filename = source["source_name"] if source else project["name"]
                    name = unique_project_name(db, filename, project["created_at"], exclude_id=project["project_id"])
                    db.execute("UPDATE projects SET name=? WHERE project_id=?", (name, project["project_id"]))
                db.execute("INSERT INTO schema_migrations(version, applied_at) VALUES (4, ?)", (utcnow(),))
        if not db.execute("SELECT 1 FROM schema_migrations WHERE version=?", (SCHEMA_VERSION,)).fetchone():
            with transaction(db):
                for project in db.execute("SELECT project_id FROM projects").fetchall():
                    project_id = project["project_id"]
                    active = db.execute(
                        "SELECT parser_version, mmp_engine_version FROM analysis_revisions WHERE project_id=? ORDER BY revision_id DESC LIMIT 1",
                        (project_id,),
                    ).fetchone()
                    if active is None or (active["parser_version"], active["mmp_engine_version"]) == (PARSER_VERSION, MMP_ENGINE_VERSION):
                        continue
                    previous = db.execute(
                        """SELECT * FROM imports WHERE project_id=? AND import_status='complete'
                           ORDER BY dataset_revision DESC LIMIT 1""",
                        (project_id,),
                    ).fetchone()
                    if previous is None or not Path(previous["source_path"]).is_file():
                        continue
                    revision = db.execute(
                        "SELECT COALESCE(MAX(dataset_revision), 0) + 1 FROM imports WHERE project_id=?",
                        (project_id,),
                    ).fetchone()[0]
                    now = utcnow()
                    import_id = db.execute(
                        """INSERT INTO imports(project_id, source_name, source_format, file_sha256,
                           source_path, imported_at, import_status, dataset_revision, summary_json)
                           VALUES (?, ?, ?, ?, ?, ?, 'queued', ?, '{}')""",
                        (
                            project_id, previous["source_name"], previous["source_format"],
                            previous["file_sha256"], previous["source_path"], now, revision,
                        ),
                    ).lastrowid
                    db.execute(
                        """INSERT INTO jobs(project_id, import_id, job_type, status, progress, created_at, updated_at)
                           VALUES (?, ?, 'import_and_analyze', 'queued', 0, ?, ?)""",
                        (project_id, import_id, now, now),
                    )
                db.execute("INSERT INTO schema_migrations(version, applied_at) VALUES (?, ?)", (SCHEMA_VERSION, utcnow()))


def row_dict(row: sqlite3.Row | None) -> dict[str, Any] | None:
    return dict(row) if row is not None else None

