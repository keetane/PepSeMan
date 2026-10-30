from __future__ import annotations

from fastapi.testclient import TestClient

from app.config import Settings
from app.chemistry import decompose_smiles

from .test_service import D3, L5, csv_bytes


def test_mmp_display_sort_handles_computed_columns_and_nulls() -> None:
    from app.api import sort_display_rows

    rows = [
        {"pair_key": "b", "convert": "Z", "delta_mw": None},
        {"pair_key": "a", "convert": "A", "delta_mw": 3.0},
        {"pair_key": "c", "convert": "A", "delta_mw": 1.0},
    ]
    assert [row["pair_key"] for row in sort_display_rows(rows, [("convert", False), ("delta_mw", True)])] == ["a", "c", "b"]
    assert [row["pair_key"] for row in sort_display_rows(rows, [("delta_mw", True)])] == ["a", "c", "b"]


def test_api_project_import_pair_and_cross_project_validation(tmp_path, monkeypatch) -> None:
    import app.db as db_module
    import app.service as service_module
    from app.api import app

    local_settings = Settings(tmp_path, tmp_path / "api.sqlite3", tmp_path / "assets")
    monkeypatch.setattr(db_module, "settings", local_settings)
    monkeypatch.setattr(service_module, "settings", local_settings)
    with TestClient(app) as client:
        project = client.post("/api/projects", json={"name": "API"}).json()
        other = client.post("/api/projects", json={"name": "Other"}).json()
        response = client.post(
            f"/api/projects/{project['project_id']}/imports",
            params={"source_name": "fixture.csv", "source_format": "csv"},
            content=csv_bytes([("L", L5, "100"), ("D", D3, "50")]),
        )
        assert response.status_code == 202
        import_id = response.json()["import_id"]
        with db_module.connect() as db:
            service_module.process_import(db, import_id)
        filtered_qc = client.get(f"/api/projects/{project['project_id']}/qc", params={"id_filter": "l", "page_size": 5})
        assert filtered_qc.status_code == 200
        assert filtered_qc.json()["total"] == 1
        assert filtered_qc.json()["rows"][0]["source_molecule_id"] == "L"
        assert client.get(f"/api/projects/{project['project_id']}/qc", params={"id_filter": "%"}).json()["total"] == 0
        assert client.get(f"/api/projects/{project['project_id']}/qc", params={"page": 2, "page_size": 1}).json()["rows"][0]["source_molecule_id"] == "D"
        decompositions = client.get(f"/api/projects/{project['project_id']}/decompositions")
        assert decompositions.status_code == 200
        assert all("value" in row and "n_pairs" in row for row in decompositions.json()["rows"])
        assert all("sequence" in row and row["sequence"] for row in decompositions.json()["rows"])
        assert max(row["n_pairs"] for row in decompositions.json()["rows"]) == 1
        unpaginated = client.get(f"/api/projects/{project['project_id']}/decompositions", params={"all_rows": "true", "page_size": 1})
        assert len(unpaginated.json()["rows"]) == 2
        filtered = client.get(f"/api/projects/{project['project_id']}/decompositions", params={"all_rows": "true", "sequence_filter": "AaA"})
        assert filtered.json()["total"] == 1
        mw_values = client.get(f"/api/projects/{project['project_id']}/decompositions", params={"value_metric": "molecular_weight"}).json()["rows"]
        assert all(row["value"] == row["molecular_weight"] for row in mw_values)
        assert client.get(f"/api/projects/{project['project_id']}/decompositions", params={"sort": "value:desc"}).status_code == 200
        pairs = client.get(f"/api/projects/{project['project_id']}/pairs", params={"sort": "delta_mw:desc"})
        assert pairs.status_code == 200
        assert pairs.json()["total"] == 1
        pair_key = pairs.json()["rows"][0]["pair_key"]
        first_hash = pairs.json()["rows"][0]["first_structure_hash"]
        second_hash = pairs.json()["rows"][0]["second_structure_hash"]
        selectable = client.get(f"/api/projects/{project['project_id']}/pairs/selection", params={"reference_hash": first_hash})
        assert selectable.status_code == 200 and selectable.json()[0]["pair_key"] == pair_key
        pair_row = pairs.json()["rows"][0]
        assert pair_row["parent_sequence"] == next(row["sequence"] for row in decompositions.json()["rows"] if row["structure_hash"] == first_hash)
        assert pair_row["parent_sequence"][pair_row["parent_sequence_highlight_start"]:pair_row["parent_sequence_highlight_end"]] == pair_row["original"]
        assert pair_row["original"] and pair_row["convert"]
        for column in ("parent_sequence", "original", "convert"):
            assert client.get(f"/api/projects/{project['project_id']}/pairs", params={"reference_hash": first_hash, "sort": f"{column}:asc"}).status_code == 200
        reversed_pairs = client.get(f"/api/projects/{project['project_id']}/pairs", params={"reference_hash": second_hash})
        assert reversed_pairs.json()["total"] == 1
        assert reversed_pairs.json()["rows"][0]["before_smiles"] == pairs.json()["rows"][0]["after_smiles"]
        assert reversed_pairs.json()["rows"][0]["delta_mw"] == -pairs.json()["rows"][0]["delta_mw"]
        assert reversed_pairs.json()["rows"][0]["original"] == pair_row["convert"]
        assert reversed_pairs.json()["rows"][0]["convert"] == pair_row["original"]
        assert client.get(f"/api/projects/{project['project_id']}/pairs", params={"reference_hash": "missing"}).status_code == 404
        detail = client.get(f"/api/projects/{project['project_id']}/pairs/{pair_key}")
        assert detail.status_code == 200
        assert "<svg" in detail.json()["before"]["svg"]
        assert "<svg" in detail.json()["scaffold_svg"]
        assert "#FFE033" in detail.json()["scaffold_svg"]
        assert detail.json()["parent_sequence"] == pair_row["parent_sequence"]
        assert detail.json()["parent_sequence"][detail.json()["parent_sequence_highlight_start"]:detail.json()["parent_sequence_highlight_end"]] == pair_row["original"]
        assert "<svg" in detail.json()["before"]["bb_svg"]
        assert "<svg" in detail.json()["after"]["bb_svg"]
        assert detail.json()["properties"]["delta_mw"] == pairs.json()["rows"][0]["delta_mw"]
        reversed_detail = client.get(f"/api/projects/{project['project_id']}/pairs/{pair_key}", params={"reference_hash": second_hash})
        assert reversed_detail.json()["before"]["structure_hash"] == second_hash
        assert reversed_detail.json()["after"]["structure_hash"] == first_hash
        assert reversed_detail.json()["properties"]["delta_mw"] == -detail.json()["properties"]["delta_mw"]
        reversed_pair_row = reversed_pairs.json()["rows"][0]
        assert reversed_detail.json()["parent_sequence"] == reversed_pair_row["parent_sequence"]
        assert reversed_pair_row["parent_sequence"][reversed_pair_row["parent_sequence_highlight_start"]:reversed_pair_row["parent_sequence_highlight_end"]] == reversed_pair_row["original"]
        selected_export = client.post(
            f"/api/projects/{project['project_id']}/export",
            json={"mode": "selected_compounds", "structure_hashes": [first_hash, second_hash]},
        )
        assert selected_export.status_code == 200
        assert first_hash in selected_export.text and second_hash in selected_export.text
        assert "sequence" in selected_export.text.splitlines()[0]
        assert client.post(
            f"/api/projects/{project['project_id']}/export",
            json={"mode": "selected_compounds", "structure_hashes": []},
        ).status_code == 422
        assert client.get(f"/api/projects/{other['project_id']}/pairs/{pair_key}").status_code == 404
        assert client.get(f"/api/projects/{project['project_id']}/pairs", params={"sort": "not_a_column:asc"}).status_code == 422


def test_sequence_dictionary_is_shared_and_editable_across_projects(tmp_path, monkeypatch) -> None:
    import app.db as db_module
    import app.service as service_module
    from app.api import app

    local_settings = Settings(tmp_path, tmp_path / "dictionary.sqlite3", tmp_path / "assets")
    monkeypatch.setattr(db_module, "settings", local_settings)
    monkeypatch.setattr(service_module, "settings", local_settings)
    n_methyl_alanine = "N[C@@H](C)C(=O)N[C@@H](C)C(=O)N(C)[C@@H](C)C(=O)N[C@@H](C)C(=O)N[C@@H](C)C(=O)O"
    with TestClient(app) as client:
        project_ids = []
        for name in ("First", "Second"):
            project_id = client.post("/api/projects", json={"name": name}).json()["project_id"]
            project_ids.append(project_id)
            import_id = client.post(
                f"/api/projects/{project_id}/imports",
                params={"source_name": f"{name}.csv", "source_format": "csv"},
                content=csv_bytes([(name, L5, "100"), (f"{name}-NMe", n_methyl_alanine, "50")]),
            ).json()["import_id"]
            with db_module.connect() as db:
                service_module.process_import(db, import_id)
        entries = client.get("/api/sequence-dictionary").json()
        alanine_hash = decompose_smiles(L5).residues[0].variant_hash
        assert all(item["variant_hash"] != alanine_hash for item in entries)
        assert all(item["default_abbreviation"].startswith("[") for item in entries)
        assert all(item["default_abbreviation"] != "[Ac]" for item in entries)
        assert client.put(f"/api/sequence-dictionary/{alanine_hash}", json={"abbreviation": "Ala"}).status_code == 404
        with db_module.connect() as db:
            db.execute("INSERT INTO sequence_aliases(variant_hash, abbreviation, updated_at) VALUES (?, 'LegacyAla', '2026-01-01')", (alanine_hash,))
        assert any(row["sequence"] == "AAAAA" for row in client.get(f"/api/projects/{project_ids[0]}/decompositions").json()["rows"])
        nonnatural = next(item for item in entries if item["variant_hash"] == decompose_smiles(n_methyl_alanine).residues[2].variant_hash)
        variant_hash = nonnatural["variant_hash"]
        svg = client.get(nonnatural["structure_url"])
        assert svg.status_code == 200 and "<svg" in svg.text
        assert client.put(f"/api/sequence-dictionary/{variant_hash}", json={"abbreviation": "NMeA"}).status_code == 200
        for project_id in project_ids:
            rows = client.get(f"/api/projects/{project_id}/decompositions").json()["rows"]
            assert any(row["sequence"] == "AA[NMeA]AA" for row in rows)
        first_rows = client.get(f"/api/projects/{project_ids[0]}/decompositions").json()["rows"]
        reference_hash = next(row["structure_hash"] for row in first_rows if row["sequence"] == "AA[NMeA]AA")
        selected_export = client.post(f"/api/projects/{project_ids[0]}/export", json={"mode": "selected_compounds", "structure_hashes": [reference_hash]})
        assert "[NMeA]" in selected_export.text
        assert client.put(f"/api/sequence-dictionary/{variant_hash}", json={"abbreviation": "bad code"}).status_code == 422
        assert client.put(f"/api/sequence-dictionary/{variant_hash}", json={"abbreviation": ""}).status_code == 200
        assert any(row["sequence"] == f"AA{nonnatural['default_abbreviation']}AA" for row in client.get(f"/api/projects/{project_ids[0]}/decompositions").json()["rows"])


def test_project_names_use_creation_date_and_source_filename(tmp_path, monkeypatch) -> None:
    import app.db as db_module
    from app.api import app

    local_settings = Settings(tmp_path, tmp_path / "names.sqlite3", tmp_path / "assets")
    monkeypatch.setattr(db_module, "settings", local_settings)
    with TestClient(app) as client:
        first = client.post("/api/projects", json={"name": "folder\\MDM2.csv"}).json()
        second = client.post("/api/projects", json={"name": "MDM2.csv"}).json()
        assert first["name"].endswith("_MDM2")
        assert second["name"] == first["name"] + "_2"
        assert len(first["name"].split("_", 1)[0]) == 8
        dated_source = client.post("/api/projects", json={"name": "20200101_experiment.csv"}).json()
        assert dated_source["name"].endswith("_20200101_experiment")
        renamed = client.patch(f"/api/projects/{second['project_id']}", json={"name": "renamed.sdf"}).json()
        assert renamed["name"].startswith(first["name"][:8] + "_renamed")


def test_existing_project_name_migration_uses_first_import_filename(tmp_path) -> None:
    from app.db import connect, initialize

    path = tmp_path / "migration.sqlite3"
    initialize(path)
    with connect(path) as db:
        db.execute("INSERT INTO projects(name, created_at, updated_at) VALUES (?, ?, ?)", ("ChEMBL MDM2 prototype", "2026-10-28T11:10:54+00:00", "2026-10-28T11:10:54+00:00"))
        db.execute("""INSERT INTO imports(project_id, source_name, source_format, file_sha256, source_path,
                      imported_at, import_status, dataset_revision) VALUES (1, 'MDM2.csv', 'csv', 'hash', 'MDM2.csv',
                      '2026-10-28T11:11:00+00:00', 'complete', 1)""")
        db.execute("DELETE FROM schema_migrations WHERE version=4")
    initialize(path)
    initialize(path)
    with connect(path) as db:
        assert db.execute("SELECT name FROM projects WHERE project_id=1").fetchone()["name"] == "20261028_MDM2"

