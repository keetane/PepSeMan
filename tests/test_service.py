from __future__ import annotations

import sqlite3
from pathlib import Path

from app.db import connect, initialize, utcnow
from app.config import Settings
from app.service import process_import, stage_import


L5 = "N[C@@H](C)C(=O)N[C@@H](C)C(=O)N[C@@H](C)C(=O)N[C@@H](C)C(=O)N[C@@H](C)C(=O)O"
D3 = "N[C@@H](C)C(=O)N[C@@H](C)C(=O)N[C@H](C)C(=O)N[C@@H](C)C(=O)N[C@@H](C)C(=O)O"
D34 = "N[C@@H](C)C(=O)N[C@@H](C)C(=O)N[C@H](C)C(=O)N[C@H](C)C(=O)N[C@@H](C)C(=O)O"


def make_db(tmp_path: Path) -> sqlite3.Connection:
    path = tmp_path / "test.sqlite3"
    initialize(path)
    return connect(path)


def project(db: sqlite3.Connection, name: str = "Test") -> int:
    now = utcnow()
    return db.execute(
        "INSERT INTO projects(name, settings_json, created_at, updated_at) VALUES (?, '{}', ?, ?)",
        (name, now, now),
    ).lastrowid


def csv_bytes(rows: list[tuple[str, str, str]]) -> bytes:
    lines = [
        '"Molecule ChEMBL ID";"Smiles";"Standard Type";"Standard Relation";"Standard Value";"Standard Units";"pChEMBL Value";"Assay ChEMBL ID"'
    ]
    for identifier, smiles, value in rows:
        lines.append(f'"{identifier}";"{smiles}";"IC50";"\'=\'";"{value}";"nM";"7.0";"ASSAY-1"')
    return ("\n".join(lines) + "\n").encode()


def test_import_keeps_measurements_separate_and_builds_strict_pairs(tmp_path: Path, monkeypatch) -> None:
    from app import service

    monkeypatch.setattr(service, "settings", Settings(tmp_path, tmp_path / "service.sqlite3", tmp_path / "assets"))
    db = make_db(tmp_path)
    project_id = project(db)
    staged = stage_import(
        db,
        project_id,
        "fixture.csv",
        "csv",
        csv_bytes([("L", L5, "100"), ("L", L5, "120"), ("D3", D3, "50"), ("D34", D34, "25")]),
    )
    summary = process_import(db, int(staged["import_id"]))
    assert summary["records"] == 4
    assert summary["unique_structures"] == 3
    assert db.execute("SELECT count(*) FROM activities WHERE project_id=?", (project_id,)).fetchone()[0] == 4
    pairs = db.execute("SELECT * FROM pairs WHERE project_id=?", (project_id,)).fetchall()
    # L5 <-> D3 and D3 <-> D34 differ at one residue; L5 <-> D34 differs at two.
    assert len(pairs) == 2
    assert {row["site_position"] for row in pairs} == {3, 4}


def test_identical_reimport_does_not_enqueue_reanalysis(tmp_path: Path, monkeypatch) -> None:
    from app import service

    monkeypatch.setattr(service, "settings", Settings(tmp_path, tmp_path / "service.sqlite3", tmp_path / "assets"))
    db = make_db(tmp_path)
    project_id = project(db)
    content = csv_bytes([("L", L5, "100")])
    first = stage_import(db, project_id, "one.csv", "csv", content)
    process_import(db, int(first["import_id"]))
    second = stage_import(db, project_id, "one-again.csv", "csv", content)
    assert second["duplicate"] is True
    assert second["job_id"] is None


def test_projects_do_not_share_pairs(tmp_path: Path, monkeypatch) -> None:
    from app import service

    monkeypatch.setattr(service, "settings", Settings(tmp_path, tmp_path / "service.sqlite3", tmp_path / "assets"))
    db = make_db(tmp_path)
    first_project = project(db, "First")
    second_project = project(db, "Second")
    staged = stage_import(db, first_project, "pair.csv", "csv", csv_bytes([("L", L5, "1"), ("D", D3, "2")]))
    process_import(db, int(staged["import_id"]))
    assert db.execute("SELECT count(*) FROM pairs WHERE project_id=?", (first_project,)).fetchone()[0] == 1
    assert db.execute("SELECT count(*) FROM pairs WHERE project_id=?", (second_project,)).fetchone()[0] == 0


def test_activity_only_update_reuses_structure_cache_and_activates_new_snapshot(tmp_path: Path, monkeypatch) -> None:
    from app import service

    monkeypatch.setattr(service, "settings", Settings(tmp_path, tmp_path / "service.sqlite3", tmp_path / "assets"))
    db = make_db(tmp_path)
    project_id = project(db)
    first = stage_import(db, project_id, "first.csv", "csv", csv_bytes([("L", L5, "100"), ("D", D3, "50")]))
    process_import(db, int(first["import_id"]))
    second = stage_import(db, project_id, "changed.csv", "csv", csv_bytes([("L", L5, "90"), ("D", D3, "40")]))
    summary = process_import(db, int(second["import_id"]))
    assert summary["recomputed"] == 0
    assert summary["affected_context_buckets"] == 0
    assert db.execute("SELECT count(*) FROM activities WHERE project_id=? AND active=1", (project_id,)).fetchone()[0] == 2
    assert db.execute("SELECT count(*) FROM activities WHERE project_id=? AND active=0", (project_id,)).fetchone()[0] == 2
    assert db.execute("SELECT count(*) FROM pairs WHERE project_id=?", (project_id,)).fetchone()[0] == 1


