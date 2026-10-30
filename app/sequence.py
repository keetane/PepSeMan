"""Conservative, deterministic peptide sequence labels from backbone decomposition."""

from __future__ import annotations

from functools import lru_cache
from typing import Mapping

from rdkit import Chem

from app.chemistry import Decomposition, Residue, decompose_smiles, sha256_text


@lru_cache(maxsize=1)
def _natural_variants() -> dict[str, str]:
    alphabet = "ACDEFGHIKLMNPQRSTVWY"
    molecule = Chem.MolFromFASTA(alphabet)
    if molecule is None:
        raise RuntimeError("RDKit could not build the natural amino acid reference")
    natural = decompose_smiles(Chem.MolToSmiles(molecule))
    if len(natural.residues) != len(alphabet):
        raise RuntimeError("Natural amino acid reference could not be decomposed")
    variants = {residue.variant_hash: letter for residue, letter in zip(natural.residues, alphabet, strict=True)}
    for residue, letter in zip(natural.residues, alphabet, strict=True):
        if letter == "G":
            continue  # Glycine has no D enantiomer.
        # Anchor indices refer to the canonicalized molecule used by the parser.
        mirror = Chem.MolFromSmiles(natural.canonical_smiles)
        # D-Ile and D-Thr also reverse the side-chain stereocenter; reversing
        # only Cα would incorrectly label the allo diastereomer.
        for atom_index in residue.atom_indices:
            atom = mirror.GetAtomWithIdx(atom_index)
            if atom.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED:
                atom.InvertChirality()
        mirrored = decompose_smiles(Chem.MolToSmiles(mirror))
        variants[mirrored.residues[residue.position - 1].variant_hash] = letter.lower()
    return variants


def _is_four_fluoroproline(residue: Residue) -> bool:
    if residue.residue_type != "Pro-like bridge" or not residue.bridge:
        return False
    bridge = Chem.MolFromSmiles(residue.bridge)
    if bridge is None:
        return False
    return sorted(atom.GetAtomicNum() for atom in bridge.GetAtoms()) == [6, 6, 6, 9]


def residue_token(value: Residue | dict, aliases: Mapping[str, str] | None = None) -> str:
    residue = Residue(**value) if isinstance(value, dict) else value
    letter = _natural_variants().get(residue.variant_hash)
    if letter is not None:
        return letter
    alias = aliases.get(residue.variant_hash) if aliases else None
    if alias:
        return f"[{alias}]"
    if _is_four_fluoroproline(residue):
        # The suffix distinguishes cis/trans and any other distinct variant.
        return f"[4FP{residue.variant_hash[:8].upper()}]"
    return f"[X{residue.variant_hash[:10].upper()}]"


def building_block_token(value: Decomposition | dict, site_kind: str, position: int, aliases: Mapping[str, str] | None = None) -> str:
    if isinstance(value, dict):
        residues = value.get("residues") or []
        n_cap, c_cap = value.get("n_cap"), value.get("c_cap")
    else:
        residues, n_cap, c_cap = value.residues, value.n_cap, value.c_cap
    if site_kind == "residue":
        return residue_token(residues[position - 1], aliases) if 0 < position <= len(residues) else ""
    if site_kind == "n_cap":
        return "[Ac]" if n_cap == "CC=O" else "[H]" if n_cap == "H" else f"[N{sha256_text(str(n_cap))[:10].upper()}]"
    if site_kind == "c_cap":
        return "[OH]" if c_cap == "O" else "[NH2]" if c_cap == "N" else f"[C{sha256_text(str(c_cap))[:10].upper()}]"
    return ""


def sequence_from_decomposition(value: Decomposition | dict, aliases: Mapping[str, str] | None = None, mask_position: int | None = None) -> str:
    if isinstance(value, dict):
        if value.get("status") not in {"accepted", "excluded_short_peptide"}:
            return ""
        n_cap = value.get("n_cap")
        residues = value.get("residues") or []
    else:
        if value.status not in {"accepted", "excluded_short_peptide"}:
            return ""
        n_cap = value.n_cap
        residues = value.residues

    parts: list[str] = []
    if n_cap == "CC=O":
        parts.append("[Ac]")
    elif n_cap not in {None, "H"}:
        parts.append(f"[N{sha256_text(n_cap)[:10].upper()}]")

    for position, raw in enumerate(residues, start=1):
        parts.append("-X-" if position == mask_position else residue_token(raw, aliases))
    return "".join(parts)


def sequence_with_site(
    value: Decomposition | dict,
    site_kind: str,
    position: int,
    aliases: Mapping[str, str] | None = None,
) -> tuple[str, int, int]:
    """Return the parent sequence and the exact span of the changed building block."""
    sequence = sequence_from_decomposition(value, aliases)
    residues = (value.get("residues") or []) if isinstance(value, dict) else value.residues
    if site_kind == "residue" and 0 < position <= len(residues):
        cap_length = len(sequence) - sum(len(residue_token(raw, aliases)) for raw in residues)
        start = cap_length + sum(len(residue_token(raw, aliases)) for raw in residues[:position - 1])
        return sequence, start, start + len(residue_token(residues[position - 1], aliases))
    if site_kind == "n_cap":
        cap = building_block_token(value, site_kind, position, aliases)
        if not sequence.startswith(cap):
            sequence = cap + sequence  # Explicitly show an otherwise implicit H cap.
        return sequence, 0, len(cap)
    if site_kind == "c_cap":
        cap = building_block_token(value, site_kind, position, aliases)
        return sequence + cap, len(sequence), len(sequence) + len(cap)
    return sequence, 0, 0
