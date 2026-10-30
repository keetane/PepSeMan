from __future__ import annotations

from pathlib import Path

import pytest
from rdkit import Chem

from app.chemistry import decompose_smiles
from app.importers import csv_to_sdf, read_csv, semantic_roundtrip_report


SOURCE_DIR = Path(__file__).parents[1] / "data" / "source"
DATASET = SOURCE_DIR / "20261030_MDM2.csv"
GLP1_DATASET = SOURCE_DIR / "20261030_GLP-1.csv"
RAP1_DATASET = SOURCE_DIR / "RAP1.csv"


@pytest.mark.skipif(not DATASET.is_file(), reason="local MDM2 source CSV is not installed")
def test_attached_chembl_csv_baseline() -> None:
    records = list(read_csv(DATASET))
    assert len(records) == 492
    assert len(records[0].properties) == 49
    assert len({item.source_molecule_id for item in records}) == 404
    assert sum(Chem.MolFromSmiles(item.smiles) is not None for item in records) == 492
    assert {item.properties["Target ChEMBL ID"] for item in records} == {"CHEMBL5023"}


@pytest.mark.skipif(not GLP1_DATASET.is_file(), reason="local GLP-1 source CSV is not installed")
def test_glp1_csv_baseline_and_compound_ids() -> None:
    records = list(read_csv(GLP1_DATASET))
    assert len(records) == 2131
    assert len({item.source_molecule_id for item in records}) == 2131
    assert all(item.source_molecule_id == item.properties["Compound ChEMBL ID"] for item in records)
    assert all(item.smiles for item in records)


@pytest.mark.skipif(not DATASET.is_file(), reason="local MDM2 source CSV is not installed")
def test_csv_sdf_roundtrip_preserves_structures_and_properties(tmp_path: Path) -> None:
    sdf = tmp_path / "mdm2.sdf"
    assert csv_to_sdf(DATASET, sdf) == 492
    report = semantic_roundtrip_report(DATASET, sdf)
    assert report["csv_records"] == report["sdf_records"] == 492
    assert report["errors"] == []


@pytest.mark.skipif(not RAP1_DATASET.is_file(), reason="local RAP1 source CSV is not installed")
def test_local_rap1_triazole_classification() -> None:
    records = list(read_csv(RAP1_DATASET))
    assert len(records) == 25
    decompositions = [decompose_smiles(record.smiles) for record in records]
    assert all(item.status == "accepted" for item in decompositions)
    assert sum(any(link["kind"] == "triazole" for link in item.cyclization_links) for item in decompositions) == 12
    assert sum(not item.cyclization_links for item in decompositions) == 13


@pytest.mark.skipif(not RAP1_DATASET.is_file(), reason="local RAP1 source CSV is not installed")
def test_local_rap1_mmp_topology(tmp_path: Path, monkeypatch) -> None:
    import app.service as service_module
    from app.config import PARSER_VERSION, Settings
    from app.db import connect, initialize
    from tests.test_service import project

    local_settings = Settings(tmp_path, tmp_path / "rap1.sqlite3", tmp_path / "assets")
    monkeypatch.setattr(service_module, "settings", local_settings)
    initialize(local_settings.database_path)
    with connect(local_settings.database_path) as db:
        project_id = project(db, "RAP1 local test")
        staged = service_module.stage_import(db, project_id, RAP1_DATASET.name, "csv", RAP1_DATASET.read_bytes())
        summary = service_module.process_import(db, int(staged["import_id"]))
        assert summary["accepted_structures"] == 25
        assert db.execute("SELECT count(*) FROM cyclization_links WHERE parser_version=? AND kind='triazole'", (PARSER_VERSION,)).fetchone()[0] == 12
        pair_counts = {
            row["kind"]: row["count"] for row in db.execute(
                """SELECT COALESCE(c.kind, 'linear') AS kind, count(*) AS count
                   FROM pairs p LEFT JOIN cyclization_links c
                     ON c.structure_hash=p.first_structure_hash AND c.parser_version=?
                   WHERE p.project_id=? GROUP BY kind""",
                (PARSER_VERSION, project_id),
            )
        }
        assert pair_counts == {"linear": 6, "triazole": 16}

