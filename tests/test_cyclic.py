"""Small synthetic structures exercise ring assignment without claiming real-data accuracy."""

from __future__ import annotations

from pathlib import Path

from rdkit import Chem
from rdkit.Chem import rdCoordGen
from fastapi.testclient import TestClient

from app.chemistry import (
    _cyclic_scaffold_molecule,
    decompose_smiles,
    parent_scaffold_svg,
    residue_building_block_svg,
)
from app.config import Settings
from app.db import SCHEMA_VERSION, connect, initialize, utcnow
from app.mmp import _matching_sites
from app.sequence import sequence_from_decomposition
from app.service import process_import, stage_import
from tests.test_service import csv_bytes, project


def peptide(sequence: str, *, head_tail: bool = False, bridges: tuple[tuple[str, int, int], ...] = ()) -> str:
    plain = Chem.MolToSmiles(Chem.MolFromSequence(sequence))
    parsed = decompose_smiles(plain)
    molecule = Chem.RWMol(Chem.MolFromSmiles(parsed.canonical_smiles))
    for kind, first, second in bridges:
        left = parsed.residues[first - 1]
        right = parsed.residues[second - 1]
        left_atoms = [index for index in left.atom_indices if molecule.GetAtomWithIdx(index).GetAtomicNum() == 16]
        if kind == "disulfide":
            right_atoms = [index for index in right.atom_indices if molecule.GetAtomWithIdx(index).GetAtomicNum() == 16]
            a, b, order = left_atoms[0], right_atoms[0], Chem.BondType.SINGLE
        elif kind == "thioether":
            # The distal Lys carbon is already linked to its terminal amine.
            right_atoms = [index for index in right.atom_indices if molecule.GetAtomWithIdx(index).GetAtomicNum() == 6
                           and any(neighbor.GetAtomicNum() == 7 for neighbor in molecule.GetAtomWithIdx(index).GetNeighbors())
                           and index != right.anchors["ca"]]
            a, b, order = left_atoms[0], right_atoms[0], Chem.BondType.SINGLE
        elif kind == "olefin_metathesis":
            def beta(residue):
                ca = residue.anchors["ca"]
                return next(atom.GetIdx() for atom in molecule.GetAtomWithIdx(ca).GetNeighbors()
                            if atom.GetAtomicNum() == 6 and atom.GetIdx() != residue.anchors["carbonyl"])
            a, b, order = beta(left), beta(right), Chem.BondType.DOUBLE
        else:
            raise ValueError(kind)
        molecule.AddBond(a, b, order)
    if head_tail:
        terminal_carbon = parsed.residues[-1].anchors["carbonyl"]
        terminal_oxygen = next(
            atom.GetIdx() for atom in molecule.GetAtomWithIdx(terminal_carbon).GetNeighbors()
            if atom.GetAtomicNum() == 8 and molecule.GetBondBetweenAtoms(terminal_carbon, atom.GetIdx()).GetBondType() == Chem.BondType.SINGLE
        )
        molecule.AddBond(terminal_carbon, parsed.residues[0].anchors["n"], Chem.BondType.SINGLE)
        molecule.RemoveAtom(terminal_oxygen)
    Chem.SanitizeMol(molecule)
    return Chem.MolToSmiles(molecule)


def triazole_peptide(sequence: str, *, close_bridge: bool = True) -> str:
    """Synthetic 1,2,3-triazole linker between two Lys side chains."""
    assert sequence[0] == sequence[-1] == "K"
    parsed = decompose_smiles(Chem.MolToSmiles(Chem.MolFromSequence(sequence)))
    molecule = Chem.RWMol(Chem.MolFromSmiles(parsed.canonical_smiles))
    terminals = [
        next(index for index in parsed.residues[position].atom_indices
             if molecule.GetAtomWithIdx(index).GetAtomicNum() == 7
             and index != parsed.residues[position].anchors["n"])
        for position in (0, -1)
    ]
    for index in terminals:
        molecule.GetAtomWithIdx(index).SetAtomicNum(6)
    base = molecule.GetMol()
    ring = Chem.MolFromSmiles("Cn1cc(C)nn1")
    match = ring.GetSubstructMatch(Chem.MolFromSmarts("[n]1[c][c][n][n]1"))
    n_arm = next(atom.GetIdx() for atom in ring.GetAtomWithIdx(match[0]).GetNeighbors() if atom.GetIdx() not in match)
    c_arm = next(atom.GetIdx() for atom in ring.GetAtomWithIdx(match[2]).GetNeighbors() if atom.GetIdx() not in match)
    combined = Chem.RWMol(Chem.CombineMols(base, ring))
    combined.AddBond(terminals[0], base.GetNumAtoms() + n_arm, Chem.BondType.SINGLE)
    if close_bridge:
        combined.AddBond(terminals[1], base.GetNumAtoms() + c_arm, Chem.BondType.SINGLE)
    else:
        # A pendant triazole is not a cyclization bond.
        combined.RemoveAtom(base.GetNumAtoms() + c_arm)
    Chem.SanitizeMol(combined)
    return Chem.MolToSmiles(combined)


def sequence_row(smiles: str) -> dict:
    result = decompose_smiles(smiles)
    assert result.status == "accepted", result.error
    return {
        "n_cap": result.n_cap,
        "c_cap": result.c_cap,
        "backbone_kind": result.backbone_kind,
        "topology_signature": result.topology_signature,
        "links": result.cyclization_links,
        "residues": [residue.variant_hash for residue in result.residues],
    }


def test_cyclic_scaffolds_use_coordgen_without_backbone_distortion() -> None:
    cases = (
        peptide("ACDEFG", head_tail=True),
        peptide("CAAAK", bridges=(("thioether", 1, 5),)),
        peptide("CAAAC", bridges=(("disulfide", 1, 5),)),
        peptide("AAAAA", bridges=(("olefin_metathesis", 1, 5),)),
    )
    for smiles in cases:
        result = decompose_smiles(smiles)
        assert result.status == "accepted"
        molecule = _cyclic_scaffold_molecule(result.canonical_smiles, result.residues)
        plain = Chem.MolFromSmiles(result.canonical_smiles)
        rdCoordGen.AddCoords(plain)
        conformer = molecule.GetConformer()
        assert all(
            (conformer.GetAtomPosition(atom.GetIdx()) - plain.GetConformer().GetAtomPosition(atom.GetIdx())).Length() < 1e-6
            for atom in molecule.GetAtoms()
        )
        assert all(
            0.8 < (conformer.GetAtomPosition(bond.GetBeginAtomIdx()) - conformer.GetAtomPosition(bond.GetEndAtomIdx())).Length() < 1.35
            for bond in molecule.GetBonds()
        )
        svg = parent_scaffold_svg(
            result.canonical_smiles, result.residues,
            changed_position=2, cyclization_links=result.cyclization_links,
        )
        assert "#FFE033" in svg


def test_head_to_tail_normalizes_sequence_and_smiles_order() -> None:
    smiles = peptide("AAVAA", head_tail=True)
    result = decompose_smiles(smiles)
    reordered = Chem.MolToSmiles(Chem.RenumberAtoms(Chem.MolFromSmiles(smiles), list(reversed(range(Chem.MolFromSmiles(smiles).GetNumAtoms())))), canonical=False)
    other = decompose_smiles(reordered)
    assert result.status == other.status == "accepted"
    assert result.backbone_kind == "head_to_tail"
    assert [link["kind"] for link in result.cyclization_links] == ["head_to_tail"]
    assert result.topology_signature == other.topology_signature
    assert sequence_from_decomposition(result) == sequence_from_decomposition(other)
    assert _matching_sites(sequence_row(peptide("AAAAA", head_tail=True)), sequence_row(smiles)) == []


def test_cyclic_rotation_and_unique_one_site_pair() -> None:
    parent = peptide("ACDEFG", head_tail=True)
    rotated = peptide("DEFGAC", head_tail=True)
    partner = peptide("ACDVFG", head_tail=True)
    assert decompose_smiles(parent).structure_hash == decompose_smiles(rotated).structure_hash
    assert len(_matching_sites(sequence_row(parent), sequence_row(partner))) == 1


def test_sidechain_bridges_and_mixed_topology() -> None:
    disulfide = decompose_smiles(peptide("CAAAC", bridges=(("disulfide", 1, 5),)))
    thioether = decompose_smiles(peptide("CAAAK", bridges=(("thioether", 1, 5),)))
    mixed = decompose_smiles(peptide("CAAACK", head_tail=True, bridges=(("disulfide", 1, 5),)))
    assert disulfide.status == thioether.status == mixed.status == "accepted"
    assert {link["kind"] for link in disulfide.cyclization_links} == {"disulfide"}
    assert {link["kind"] for link in thioether.cyclization_links} == {"thioether"}
    assert {link["kind"] for link in mixed.cyclization_links} == {"head_to_tail", "disulfide"}
    assert _matching_sites(sequence_row(peptide("CAAAC")), sequence_row(peptide("CAAAC", bridges=(("disulfide", 1, 5),)))) == []


def test_ring_alkene_is_structural_category_not_plain_alkene() -> None:
    ring = decompose_smiles(peptide("AAAAA", bridges=(("olefin_metathesis", 1, 5),)))
    plain = decompose_smiles(Chem.MolToSmiles(Chem.MolFromSequence("AAAAA")))
    assert ring.status == plain.status == "accepted"
    assert {link["kind"] for link in ring.cyclization_links} == {"olefin_metathesis"}
    assert plain.cyclization_links == []
    assert ring.topology_signature != plain.topology_signature
    assert not decompose_smiles(Chem.MolToSmiles(Chem.MolFromSequence("AAMAA"))).cyclization_links


def test_triazole_bridge_is_structural_and_mmp_topology_is_strict() -> None:
    parent = triazole_peptide("KAAAK")
    partner = triazole_peptide("KAVAK")
    open_sidechain = triazole_peptide("KAAAK", close_bridge=False)
    a, b, pendant = (decompose_smiles(smiles) for smiles in (parent, partner, open_sidechain))
    assert a.status == b.status == pendant.status == "accepted"
    assert [link["kind"] for link in a.cyclization_links] == ["triazole"]
    assert a.cyclization_links[0]["linker_smiles"].count("n") == 3
    assert a.topology_signature == b.topology_signature
    assert len(_matching_sites(sequence_row(parent), sequence_row(partner))) == 1
    assert pendant.cyclization_links == []
    assert pendant.topology_signature != a.topology_signature
    assert _matching_sites(sequence_row(parent), sequence_row(open_sidechain)) == []


def test_triazole_smiles_order_does_not_change_topology() -> None:
    smiles = triazole_peptide("KAAAK")
    molecule = Chem.MolFromSmiles(smiles)
    reordered = Chem.MolToSmiles(Chem.RenumberAtoms(molecule, list(reversed(range(molecule.GetNumAtoms())))), canonical=False)
    first, second = decompose_smiles(smiles), decompose_smiles(reordered)
    assert first.status == second.status == "accepted"
    assert first.topology_signature == second.topology_signature
    assert first.cyclization_links[0]["linker_smiles"] == second.cyclization_links[0]["linker_smiles"]


def test_plain_sidechain_alkene_is_not_a_bridge() -> None:
    base = decompose_smiles(Chem.MolToSmiles(Chem.MolFromSequence("AAAAA")))
    molecule = Chem.RWMol(Chem.MolFromSmiles(base.canonical_smiles))
    residue = base.residues[2]
    ca = residue.anchors["ca"]
    beta = next(atom.GetIdx() for atom in molecule.GetAtomWithIdx(ca).GetNeighbors()
                if atom.GetAtomicNum() == 6 and atom.GetIdx() != residue.anchors["carbonyl"])
    terminal = molecule.AddAtom(Chem.Atom(6))
    molecule.AddBond(beta, terminal, Chem.BondType.DOUBLE)
    Chem.SanitizeMol(molecule)
    result = decompose_smiles(Chem.MolToSmiles(molecule))
    assert result.status == "accepted"
    assert result.cyclization_links == []


def test_non_unique_short_thioether_cut_goes_to_qc() -> None:
    base = decompose_smiles(Chem.MolToSmiles(Chem.MolFromSequence("CAAAA")))
    molecule = Chem.RWMol(Chem.MolFromSmiles(base.canonical_smiles))
    sulfur = next(index for index in base.residues[0].atom_indices if molecule.GetAtomWithIdx(index).GetAtomicNum() == 16)
    last = base.residues[-1]
    beta = next(atom.GetIdx() for atom in molecule.GetAtomWithIdx(last.anchors["ca"]).GetNeighbors()
                if atom.GetAtomicNum() == 6 and atom.GetIdx() != last.anchors["carbonyl"])
    molecule.AddBond(sulfur, beta, Chem.BondType.SINGLE)
    Chem.SanitizeMol(molecule)
    result = decompose_smiles(Chem.MolToSmiles(molecule))
    assert result.status == "unsupported_backbone"
    assert "unique" in result.error


def test_ring_alkene_e_z_is_retained_in_topology() -> None:
    base = Chem.MolFromSmiles(peptide("AAAAAA", bridges=(("olefin_metathesis", 1, 6),)))
    link = decompose_smiles(Chem.MolToSmiles(base)).cyclization_links[0]
    bond = base.GetBondBetweenAtoms(*link["bond_indices"])
    left = next(atom.GetIdx() for atom in bond.GetBeginAtom().GetNeighbors() if atom.GetIdx() != bond.GetEndAtomIdx())
    right = next(atom.GetIdx() for atom in bond.GetEndAtom().GetNeighbors() if atom.GetIdx() != bond.GetBeginAtomIdx())
    bond.SetStereoAtoms(left, right)
    signatures = []
    for stereo in (Chem.BondStereo.STEREOZ, Chem.BondStereo.STEREOE):
        bond.SetStereo(stereo)
        result = decompose_smiles(Chem.MolToSmiles(base, isomericSmiles=True))
        assert result.status == "accepted"
        assert result.cyclization_links[0]["stereo"] == str(stereo)
        signatures.append(result.topology_signature)
    assert signatures[0] != signatures[1]


def test_two_same_type_bridges_are_retained() -> None:
    result = decompose_smiles(peptide("CAAACCAAAC", bridges=(("disulfide", 1, 5), ("disulfide", 6, 10))))
    assert result.status == "accepted", result.error
    assert [link["kind"] for link in result.cyclization_links] == ["disulfide", "disulfide"]


def test_mixed_disulfide_and_ring_alkene() -> None:
    result = decompose_smiles(peptide("CAAACAAAAA", bridges=(("disulfide", 1, 5), ("olefin_metathesis", 6, 10))))
    assert result.status == "accepted", result.error
    assert {link["kind"] for link in result.cyclization_links} == {"disulfide", "olefin_metathesis"}


def test_cap_to_sidechain_thioether_has_a_cap_endpoint() -> None:
    base = decompose_smiles(Chem.MolToSmiles(Chem.MolFromSequence("AAAAC")))
    molecule = Chem.RWMol(Chem.MolFromSmiles(base.canonical_smiles))
    carbonyl = molecule.AddAtom(Chem.Atom(6))
    oxygen = molecule.AddAtom(Chem.Atom(8))
    cap_carbon = molecule.AddAtom(Chem.Atom(6))
    molecule.AddBond(base.residues[0].anchors["n"], carbonyl, Chem.BondType.SINGLE)
    molecule.AddBond(carbonyl, oxygen, Chem.BondType.DOUBLE)
    molecule.AddBond(carbonyl, cap_carbon, Chem.BondType.SINGLE)
    sulfur = next(index for index in base.residues[-1].atom_indices if molecule.GetAtomWithIdx(index).GetAtomicNum() == 16)
    molecule.AddBond(cap_carbon, sulfur, Chem.BondType.SINGLE)
    Chem.SanitizeMol(molecule)
    result = decompose_smiles(Chem.MolToSmiles(molecule))
    assert result.status == "accepted", result.error
    assert result.cyclization_links[0]["kind"] == "thioether"
    assert {result.cyclization_links[0]["first"]["kind"], result.cyclization_links[0]["second"]["kind"]} == {"n_cap", "residue"}


def test_nonnatural_bridged_building_block_marks_attachment() -> None:
    result = decompose_smiles(peptide("CAAAK", bridges=(("thioether", 1, 5),)))
    link = result.cyclization_links[0]
    attachment = next(end["atom_index"] for end in (link["first"], link["second"]) if end["position"] == 5)
    svg = residue_building_block_svg(result.canonical_smiles, result.residues[4], attachment_atoms=[attachment])
    assert "class='note'" in svg


def test_upgrade_queues_reanalysis_without_replacing_active_revision(tmp_path: Path) -> None:
    path = tmp_path / "migration.sqlite3"
    source = tmp_path / "source.csv"
    source.write_bytes(csv_bytes([("A", peptide("CAAAC"), "100")]))
    initialize(path)
    with connect(path) as db:
        project_id = project(db)
        now = utcnow()
        db.execute(
            """INSERT INTO imports(project_id, source_name, source_format, file_sha256, source_path,
               imported_at, import_status, dataset_revision) VALUES (?, 'source.csv', 'csv', 'hash', ?, ?, 'complete', 1)""",
            (project_id, str(source), now),
        )
        db.execute(
            """INSERT INTO analysis_revisions(project_id, input_revision, parser_version,
               mmp_engine_version, summary_json, created_at) VALUES (?, 1, 'backbone-v1', 'strict-one-site-v1', '{}', ?)""",
            (project_id, now),
        )
        db.execute("DELETE FROM schema_migrations WHERE version=?", (SCHEMA_VERSION,))
    initialize(path)
    initialize(path)
    with connect(path) as db:
        assert db.execute("SELECT parser_version FROM analysis_revisions WHERE project_id=?", (project_id,)).fetchone()[0] == "backbone-v1"
        queued = db.execute("SELECT * FROM jobs WHERE project_id=?", (project_id,)).fetchall()
        assert len(queued) == 1 and queued[0]["status"] == "queued"
        import_id = queued[0]["import_id"]
        assert db.execute("SELECT dataset_revision FROM imports WHERE import_id=?", (import_id,)).fetchone()[0] == 2


def test_cyclic_project_import_preserves_topology_and_atomic_revision(tmp_path: Path, monkeypatch) -> None:
    from app import service

    monkeypatch.setattr(service, "settings", Settings(tmp_path, tmp_path / "test.sqlite3", tmp_path / "assets"))
    path = tmp_path / "test.sqlite3"
    initialize(path)
    with connect(path) as db:
        project_id = project(db)
        first = peptide("CAAAC", bridges=(("disulfide", 1, 5),))
        second = peptide("CAA VC".replace(" ", ""), bridges=(("disulfide", 1, 5),))
        staged = stage_import(db, project_id, "cyclic.csv", "csv", csv_bytes([("A", first, "100"), ("B", second, "50")]))
        summary = process_import(db, int(staged["import_id"]))
        assert summary["accepted_structures"] == 2
        assert db.execute("SELECT count(*) FROM cyclization_links WHERE kind='disulfide'").fetchone()[0] == 2
        assert db.execute("SELECT count(*) FROM pairs WHERE project_id=?", (project_id,)).fetchone()[0] == 1


def test_triazole_project_pair_is_indexed_without_open_chain_pair(tmp_path: Path, monkeypatch) -> None:
    import app.db as db_module
    from app import service
    from app.api import app

    local_settings = Settings(tmp_path, tmp_path / "test.sqlite3", tmp_path / "assets")
    monkeypatch.setattr(service, "settings", local_settings)
    monkeypatch.setattr(db_module, "settings", local_settings)
    path = tmp_path / "test.sqlite3"
    initialize(path)
    with connect(path) as db:
        project_id = project(db)
        staged = stage_import(db, project_id, "synthetic-triazole.csv", "csv", csv_bytes([
            ("A", triazole_peptide("KAAAK"), "100"),
            ("B", triazole_peptide("KAVAK"), "50"),
            ("C", triazole_peptide("KAAAK", close_bridge=False), "150"),
        ]))
        summary = process_import(db, int(staged["import_id"]))
        assert summary["accepted_structures"] == 3
        assert db.execute("SELECT count(*) FROM cyclization_links WHERE kind='triazole'").fetchone()[0] == 2
        assert db.execute("SELECT count(*) FROM pairs WHERE project_id=?", (project_id,)).fetchone()[0] == 1
    with TestClient(app) as client:
        filtered = client.get(f"/api/projects/{project_id}/decompositions", params={"cyclization_type": "triazole"})
        assert filtered.status_code == 200 and filtered.json()["total"] == 2
        reference = filtered.json()["rows"][0]["structure_hash"]
        pairs = client.get(f"/api/projects/{project_id}/pairs", params={"reference_hash": reference})
        assert pairs.status_code == 200 and pairs.json()["total"] == 1
        detail = client.get(f"/api/projects/{project_id}/pairs/{pairs.json()['rows'][0]['pair_key']}", params={"reference_hash": reference})
        assert detail.status_code == 200
        assert detail.json()["cyclization_links"][0]["kind"] == "triazole"


def test_head_to_tail_project_pair_is_found_after_context_indexing(tmp_path: Path, monkeypatch) -> None:
    from app import service

    monkeypatch.setattr(service, "settings", Settings(tmp_path, tmp_path / "test.sqlite3", tmp_path / "assets"))
    path = tmp_path / "test.sqlite3"
    initialize(path)
    with connect(path) as db:
        project_id = project(db)
        staged = stage_import(db, project_id, "rings.csv", "csv", csv_bytes([
            ("A", peptide("ACDEFG", head_tail=True), "100"),
            ("B", peptide("ACDVFG", head_tail=True), "50"),
        ]))
        process_import(db, int(staged["import_id"]))
        pair = db.execute("SELECT site_position, second_site_position FROM pairs WHERE project_id=?", (project_id,)).fetchone()
        assert pair is not None
        assert 1 <= pair["site_position"] <= 6
        assert 1 <= pair["second_site_position"] <= 6


def test_cyclic_api_filter_pair_detail_and_csv(tmp_path: Path, monkeypatch) -> None:
    import app.db as db_module
    import app.service as service_module
    from app.api import app

    local_settings = Settings(tmp_path, tmp_path / "api.sqlite3", tmp_path / "assets")
    monkeypatch.setattr(db_module, "settings", local_settings)
    monkeypatch.setattr(service_module, "settings", local_settings)
    first = peptide("CAAAC", bridges=(("disulfide", 1, 5),))
    second = peptide("CAVAC", bridges=(("disulfide", 1, 5),))
    linear = peptide("CAAAC")
    with TestClient(app) as client:
        project_id = client.post("/api/projects", json={"name": "cyclic"}).json()["project_id"]
        response = client.post(
            f"/api/projects/{project_id}/imports",
            params={"source_name": "cyclic.csv", "source_format": "csv"},
            content=csv_bytes([("A", first, "100"), ("B", second, "50"), ("L", linear, "75")]),
        )
        with db_module.connect() as db:
            service_module.process_import(db, response.json()["import_id"])
        decomposition = client.get(f"/api/projects/{project_id}/decompositions", params={"cyclization_type": "disulfide", "sort": "link_count:desc"})
        assert decomposition.status_code == 200
        assert decomposition.json()["total"] == 2
        row = decomposition.json()["rows"][0]
        assert row["cyclization_types"] == "disulfide"
        assert row["link_count"] == 1
        assert row["cyclization_links"][0]["first"]["position"] == 1
        assert client.get(f"/api/projects/{project_id}/decompositions", params={"cyclization_type": "linear"}).json()["total"] == 1
        pairs = client.get(f"/api/projects/{project_id}/pairs", params={"reference_hash": row["structure_hash"], "sort": "cyclization_types:asc"})
        assert pairs.status_code == 200 and pairs.json()["total"] == 1
        detail = client.get(f"/api/projects/{project_id}/pairs/{pairs.json()['rows'][0]['pair_key']}", params={"reference_hash": row["structure_hash"]})
        assert detail.status_code == 200
        assert "#FFE033" in detail.json()["scaffold_svg"]
        assert detail.json()["cyclization_links"][0]["kind"] == "disulfide"
        exported = client.post(f"/api/projects/{project_id}/export", json={"mode": "selected_compounds", "structure_hashes": [row["structure_hash"]]})
        assert exported.status_code == 200
        assert "cyclization_links_json" in exported.text.splitlines()[0]
        assert "disulfide" in exported.text
