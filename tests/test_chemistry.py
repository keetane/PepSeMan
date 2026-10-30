from __future__ import annotations

from rdkit import Chem

from app.chemistry import _linear_scaffold_molecule, canonicalize_smiles, decompose_smiles, parent_scaffold_svg
from app.sequence import sequence_from_decomposition, sequence_with_site


def glycine_peptide(length: int) -> str:
    return "NCC(=O)" + "NCC(=O)" * (length - 1) + "O"


def test_scaffold_backbone_runs_from_n_to_c_and_highlights_changed_block() -> None:
    smiles = glycine_peptide(5)
    residues = decompose_smiles(smiles).residues
    molecule = _linear_scaffold_molecule(smiles, residues)
    conformer = molecule.GetConformer()
    positions = [conformer.GetAtomPosition(residue.anchors["n"]).x for residue in residues]
    assert positions == sorted(positions)
    assert positions[0] < conformer.GetAtomPosition(residues[-1].anchors["carbonyl"]).x
    assert "#FFE033" in parent_scaffold_svg(smiles, residues, changed_position=3)
    assert "#FFE033" in parent_scaffold_svg(smiles, residues, site_kind="n_cap")


def test_scaffold_backbone_zigzags_and_side_chains_alternate() -> None:
    smiles = Chem.MolToSmiles(Chem.MolFromSequence("AAAAAA"))
    residues = decompose_smiles(smiles).residues
    molecule = _linear_scaffold_molecule(smiles, residues)
    conformer = molecule.GetConformer()
    side_chain_directions = []
    backbone_positions = []
    for index, residue in enumerate(residues):
        n = conformer.GetAtomPosition(residue.anchors["n"])
        ca = conformer.GetAtomPosition(residue.anchors["ca"])
        carbonyl = conformer.GetAtomPosition(residue.anchors["carbonyl"])
        assert n.x < ca.x < carbonyl.x
        backbone_positions.extend((n, ca, carbonyl))
        assert (ca.y > n.y) == (index % 2 == 0)
        side_chain_atoms = [
            neighbor.GetIdx()
            for neighbor in molecule.GetAtomWithIdx(residue.anchors["ca"]).GetNeighbors()
            if neighbor.GetIdx() not in {residue.anchors["n"], residue.anchors["carbonyl"]}
        ]
        assert len(side_chain_atoms) == 1
        side_chain = conformer.GetAtomPosition(side_chain_atoms[0])
        side_chain_directions.append(side_chain.y - ca.y)
    for index, point in enumerate(backbone_positions):
        assert abs(point.y - (0.75 if index % 2 else 0)) < 1e-6
        if index:
            previous = backbone_positions[index - 1]
            assert point.x > previous.x
            assert abs((point.x - previous.x) - 1.3) < 1e-6
    assert all(direction > 0 if index % 2 == 0 else direction < 0 for index, direction in enumerate(side_chain_directions))


def test_parent_sequence_highlights_exact_residue_or_terminal_cap() -> None:
    plain = decompose_smiles(Chem.MolToSmiles(Chem.MolFromSequence("AAPAA")))
    sequence, start, end = sequence_with_site(plain, "residue", 3)
    assert sequence == "AAPAA"
    assert sequence[start:end] == "P"
    sequence, start, end = sequence_with_site(plain, "n_cap", 1)
    assert sequence == "[H]AAPAA"
    assert sequence[start:end] == "[H]"
    sequence, start, end = sequence_with_site(plain, "c_cap", 5)
    assert sequence == "AAPAA[OH]"
    assert sequence[start:end] == "[OH]"

    fluoroproline = decompose_smiles(Chem.MolToSmiles(Chem.MolFromSequence("AAPAA")).replace("N1CCC[C@H]1", "N1C[C@H](F)C[C@H]1"))
    sequence, start, end = sequence_with_site(fluoroproline, "residue", 3)
    assert sequence[start:end].startswith("[4FP")


def test_only_five_or_more_residues_are_accepted() -> None:
    for length in range(1, 5):
        result = decompose_smiles(glycine_peptide(length))
        assert result.status == "excluded_short_peptide"
        assert result.length == length
    result = decompose_smiles(glycine_peptide(5))
    assert result.status == "accepted"
    assert result.length == 5
    assert [residue.residue_type for residue in result.residues] == ["Gly"] * 5


def test_caps_are_not_counted_and_l_stereo_is_structural() -> None:
    smiles = "CC(=O)N[C@@H](C)C(=O)N[C@@H](C)C(=O)N[C@@H](C)C(=O)N[C@@H](C)C(=O)N[C@@H](C)C(N)=O"
    result = decompose_smiles(smiles)
    assert result.status == "accepted"
    assert result.length == 5
    assert result.n_cap != "H"
    assert result.c_cap != "unspecified"
    assert [item.stereo_status for item in result.residues] == ["L"] * 5
    assert sequence_from_decomposition(result) == "[Ac]AAAAA"


def test_sequence_distinguishes_natural_d_and_non_natural_building_blocks() -> None:
    molecule = Chem.MolFromSequence("AAPAA")
    natural = decompose_smiles(Chem.MolToSmiles(molecule))
    molecule = Chem.MolFromSmiles(natural.canonical_smiles)
    molecule.GetAtomWithIdx(natural.residues[1].anchors["ca"]).InvertChirality()
    d_peptide = decompose_smiles(Chem.MolToSmiles(molecule))
    assert sequence_from_decomposition(d_peptide) == "AaPAA"

    base_smiles = Chem.MolToSmiles(Chem.MolFromSequence("AAPAA"))
    first = decompose_smiles(base_smiles.replace("N1CCC[C@H]1", "N1C[C@H](F)C[C@H]1"))
    second = decompose_smiles(base_smiles.replace("N1CCC[C@H]1", "N1C[C@@H](F)C[C@H]1"))
    first_sequence = sequence_from_decomposition(first)
    second_sequence = sequence_from_decomposition(second)
    assert first_sequence.startswith("AA[4FP") and first_sequence.endswith("]AA")
    assert second_sequence.startswith("AA[4FP") and second_sequence != first_sequence
    assert sequence_from_decomposition(first.as_json()) == first_sequence


def test_d_isoleucine_and_threonine_require_side_chain_stereo_inversion() -> None:
    for letter in "IT":
        original = decompose_smiles(Chem.MolToSmiles(Chem.MolFromSequence(f"AA{letter}AA")))
        mirrored = Chem.MolFromSmiles(original.canonical_smiles)
        for atom_index in original.residues[2].atom_indices:
            atom = mirrored.GetAtomWithIdx(atom_index)
            if atom.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED:
                atom.InvertChirality()
        assert sequence_from_decomposition(decompose_smiles(Chem.MolToSmiles(mirrored))) == f"AA{letter.lower()}AA"


def test_l_to_d_change_is_recorded_at_one_residue() -> None:
    l_peptide = "N[C@@H](C)C(=O)N[C@@H](C)C(=O)N[C@@H](C)C(=O)N[C@@H](C)C(=O)N[C@@H](C)C(=O)O"
    d_peptide = "N[C@@H](C)C(=O)N[C@@H](C)C(=O)N[C@H](C)C(=O)N[C@@H](C)C(=O)N[C@@H](C)C(=O)O"
    left = decompose_smiles(l_peptide)
    right = decompose_smiles(d_peptide)
    assert left.status == right.status == "accepted"
    differences = [
        index
        for index, (a, b) in enumerate(zip(left.residues, right.residues, strict=True), start=1)
        if a.variant_hash != b.variant_hash
    ]
    assert differences == [3]
    assert left.residues[2].stereo_status != right.residues[2].stereo_status


def test_side_chain_amide_does_not_add_a_residue() -> None:
    # Gln-Gly-Gly-Gly-Gly; the side-chain amide is not part of the backbone.
    smiles = "N[C@@H](CCC(N)=O)C(=O)NCC(=O)NCC(=O)NCC(=O)NCC(=O)O"
    result = decompose_smiles(smiles)
    assert result.status == "accepted"
    assert result.length == 5


def test_multiple_fragments_are_not_silently_desalted() -> None:
    result = decompose_smiles(glycine_peptide(5) + ".[Na+]")
    assert result.status == "unsupported_multiple_fragments"


def test_equivalent_smiles_have_identical_structure_hashes() -> None:
    _, _, first = canonicalize_smiles("NCC(=O)O")
    _, _, second = canonicalize_smiles("O=C(O)CN")
    assert first == second


def test_invalid_structure_is_reported_not_raised() -> None:
    result = decompose_smiles("this-is-not-smiles")
    assert result.status == "invalid_structure"
    assert result.structure_hash is None


def test_proline_is_an_indivisible_bridge_and_d_pro_changes_one_residue() -> None:
    molecule = Chem.MolFromSequence("AAPAA")
    l_pro = Chem.MolToSmiles(molecule, canonical=True, isomericSmiles=True)
    d_pro = l_pro.replace("N1CCC[C@H]1", "N1CCC[C@@H]1")
    left = decompose_smiles(l_pro)
    right = decompose_smiles(d_pro)
    assert left.status == right.status == "accepted"
    assert left.residues[2].residue_type == "Pro-like bridge"
    assert left.residues[2].bridge is not None
    assert left.residues[2].r1 is left.residues[2].r3 is None
    assert sum(a.variant_hash != b.variant_hash for a, b in zip(left.residues, right.residues, strict=True)) == 1
    hydroxyproline = l_pro.replace("N1CCC[C@H]1", "N1CC(O)C[C@H]1")
    hyp = decompose_smiles(hydroxyproline)
    assert hyp.status == "accepted"
    assert hyp.residues[2].residue_type == "Pro-like bridge"
    assert sum(a.variant_hash != b.variant_hash for a, b in zip(left.residues, hyp.residues, strict=True)) == 1


def test_n_methyl_peptoid_and_alpha_disubstitution_are_distinguished() -> None:
    n_methyl = "N[C@@H](C)C(=O)N[C@@H](C)C(=O)N(C)[C@@H](C)C(=O)N[C@@H](C)C(=O)N[C@@H](C)C(=O)O"
    peptoid = "N[C@@H](C)C(=O)N[C@@H](C)C(=O)N(CC)CC(=O)N[C@@H](C)C(=O)N[C@@H](C)C(=O)O"
    disubstituted = "N[C@@H](C)C(=O)N[C@@H](C)C(=O)NC(C)(C)C(=O)N[C@@H](C)C(=O)N[C@@H](C)C(=O)O"
    first = decompose_smiles(n_methyl)
    second = decompose_smiles(peptoid)
    third = decompose_smiles(disubstituted)
    assert first.residues[2].residue_type == "N-substituted amino acid"
    assert first.residues[2].r3 is not None
    assert second.residues[2].residue_type == "N-substituted Gly / peptoid"
    assert second.residues[2].r3 is not None
    assert third.residues[2].residue_type == "alpha,alpha-disubstituted"


def test_lysine_and_serine_side_chain_modifications_change_only_their_residue() -> None:
    lys = Chem.MolToSmiles(Chem.MolFromSequence("AAKAA"), canonical=True, isomericSmiles=True)
    ac_lys = lys.replace("CCCCN", "CCCCNC(C)=O")
    ser = Chem.MolToSmiles(Chem.MolFromSequence("AASAA"), canonical=True, isomericSmiles=True)
    methyl_ser = ser.replace("CO)", "COC)")
    for before, after in ((lys, ac_lys), (ser, methyl_ser)):
        left = decompose_smiles(before)
        right = decompose_smiles(after)
        assert left.status == right.status == "accepted"
        assert [
            index
            for index, (a, b) in enumerate(zip(left.residues, right.residues, strict=True), start=1)
            if a.variant_hash != b.variant_hash
        ] == [3]


def test_alanine_to_valine_is_one_residue_change() -> None:
    alanine = decompose_smiles(Chem.MolToSmiles(Chem.MolFromSequence("AAAAA"), True))
    valine = decompose_smiles(Chem.MolToSmiles(Chem.MolFromSequence("AAVAA"), True))
    assert alanine.status == valine.status == "accepted"
    assert [
        index
        for index, (a, b) in enumerate(zip(alanine.residues, valine.residues, strict=True), start=1)
        if a.variant_hash != b.variant_hash
    ] == [3]

