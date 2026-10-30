from __future__ import annotations

import hashlib
import json
from collections import deque
from dataclasses import asdict, dataclass, field
from itertools import combinations
from typing import Iterable, Sequence

from rdkit import Chem
from rdkit.Chem import Descriptors, rdCoordGen, rdDepictor, rdMolDescriptors
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Geometry import Point2D

from app.config import PARSER_VERSION, STANDARDIZATION_VERSION


# The substituted N of a 1,2,3-triazole is the first atom. Cutting its
# exocyclic N–C bond leaves the aromatic ring intact on the other residue.
_TRIAZOLE_QUERY = Chem.MolFromSmarts("[n]1[c][c][n][n]1")


@dataclass(frozen=True)
class BackboneUnit:
    n: int
    ca: int
    carbonyl: int
    oxygens: tuple[int, ...]


@dataclass
class Residue:
    position: int
    residue_type: str
    r1: str | None
    r2: str | None
    r3: str | None
    bridge: str | None
    stereo_status: str
    atom_indices: list[int]
    anchors: dict[str, int | list[int]]
    variant_hash: str


@dataclass
class Decomposition:
    status: str
    canonical_smiles: str | None
    structure_hash: str | None
    length: int | None
    n_cap: str | None
    c_cap: str | None
    residues: list[Residue]
    error: str | None = None
    molecular_weight: float | None = None
    parser_version: str = PARSER_VERSION
    standardization_version: str = STANDARDIZATION_VERSION
    backbone_kind: str = "linear"
    cyclization_links: list[dict] = field(default_factory=list)
    topology_signature: str = ""
    canonical_start: int = 1

    def as_json(self) -> dict:
        value = asdict(self)
        value["residues"] = [asdict(residue) for residue in self.residues]
        return value


def sha256_text(value: str) -> str:
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


def canonicalize_smiles(smiles: str) -> tuple[Chem.Mol, str, str]:
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        raise ValueError("RDKit could not parse SMILES")
    Chem.AssignStereochemistry(molecule, cleanIt=True, force=True)
    canonical = Chem.MolToSmiles(molecule, canonical=True, isomericSmiles=True)
    structure_hash = sha256_text(f"{STANDARDIZATION_VERSION}\0{canonical}")
    return molecule, canonical, structure_hash


def _carbonyl_oxygens(molecule: Chem.Mol, atom_index: int) -> tuple[int, ...]:
    atom = molecule.GetAtomWithIdx(atom_index)
    values: list[int] = []
    for bond in atom.GetBonds():
        other = bond.GetOtherAtom(atom)
        if other.GetAtomicNum() == 8 and bond.GetBondType() == Chem.BondType.DOUBLE:
            values.append(other.GetIdx())
    return tuple(sorted(values))


def _candidate_units(molecule: Chem.Mol) -> list[BackboneUnit]:
    candidates: set[BackboneUnit] = set()
    for carbonyl in molecule.GetAtoms():
        if carbonyl.GetAtomicNum() != 6:
            continue
        oxygens = _carbonyl_oxygens(molecule, carbonyl.GetIdx())
        if not oxygens:
            continue
        for c_bond in carbonyl.GetBonds():
            ca = c_bond.GetOtherAtom(carbonyl)
            if c_bond.GetBondType() != Chem.BondType.SINGLE or ca.GetAtomicNum() != 6:
                continue
            if ca.GetHybridization() not in (Chem.HybridizationType.SP3, Chem.HybridizationType.UNSPECIFIED):
                continue
            for n_bond in ca.GetBonds():
                nitrogen = n_bond.GetOtherAtom(ca)
                if (
                    nitrogen.GetIdx() == carbonyl.GetIdx()
                    or n_bond.GetBondType() != Chem.BondType.SINGLE
                    or nitrogen.GetAtomicNum() != 7
                ):
                    continue
                candidates.add(
                    BackboneUnit(
                        n=nitrogen.GetIdx(),
                        ca=ca.GetIdx(),
                        carbonyl=carbonyl.GetIdx(),
                        oxygens=oxygens,
                    )
                )
    return sorted(candidates, key=lambda unit: (unit.n, unit.ca, unit.carbonyl))


def _unit_graph(molecule: Chem.Mol, units: Sequence[BackboneUnit]) -> tuple[dict[int, list[int]], dict[int, list[int]]]:
    outgoing: dict[int, list[int]] = {index: [] for index in range(len(units))}
    incoming: dict[int, list[int]] = {index: [] for index in range(len(units))}
    for left_index, left in enumerate(units):
        for right_index, right in enumerate(units):
            if left_index == right_index or left.ca == right.ca:
                continue
            bond = molecule.GetBondBetweenAtoms(left.carbonyl, right.n)
            if bond is not None and bond.GetBondType() == Chem.BondType.SINGLE:
                outgoing[left_index].append(right_index)
                incoming[right_index].append(left_index)
    return outgoing, incoming


def _maximal_paths(outgoing: dict[int, list[int]], incoming: dict[int, list[int]]) -> list[tuple[int, ...]]:
    starts = [index for index, previous in incoming.items() if not previous]
    if not starts:
        return []
    results: set[tuple[int, ...]] = set()

    def walk(current: int, path: tuple[int, ...]) -> None:
        next_values = [item for item in outgoing[current] if item not in path]
        if not next_values:
            results.add(path)
            return
        for item in next_values:
            walk(item, (*path, item))

    for start in starts:
        walk(start, (start,))
    return sorted(results, key=lambda path: (-len(path), path))


def _backbone_cycles(outgoing: dict[int, list[int]], units: Sequence[BackboneUnit]) -> list[tuple[int, ...]]:
    """Find directed peptide-bond cycles, treating rotations as one candidate."""
    cycles: set[tuple[int, ...]] = set()
    for start in range(len(units)):
        stack = [(start, (start,))]
        while stack:
            current, path = stack.pop()
            for following in outgoing[current]:
                if following == start and len(path) >= 2:
                    rotations = [path[index:] + path[:index] for index in range(len(path))]
                    cycles.add(min(rotations))
                elif following not in path and len(path) < len(units):
                    stack.append((following, (*path, following)))
    return sorted(cycles, key=lambda path: (-len(path), path))


def _link_kind(bond: Chem.Bond) -> str | None:
    numbers = sorted((bond.GetBeginAtom().GetAtomicNum(), bond.GetEndAtom().GetAtomicNum()))
    if numbers == [6, 7] and bond.GetBondType() == Chem.BondType.SINGLE:
        nitrogen = bond.GetBeginAtom() if bond.GetBeginAtom().GetAtomicNum() == 7 else bond.GetEndAtom()
        carbon = bond.GetEndAtom() if nitrogen == bond.GetBeginAtom() else bond.GetBeginAtom()
        if nitrogen.GetIsAromatic() and any(
            match[0] == nitrogen.GetIdx() and carbon.GetIdx() not in match
            for match in bond.GetOwningMol().GetSubstructMatches(_TRIAZOLE_QUERY)
        ):
            return "triazole"
    if numbers == [16, 16] and bond.GetBondType() == Chem.BondType.SINGLE:
        return "disulfide"
    if numbers == [6, 16] and bond.GetBondType() == Chem.BondType.SINGLE:
        # A C-S bond of cystine is part of a disulfide side chain, not a
        # competing thioether cut.
        sulfur = bond.GetBeginAtom() if bond.GetBeginAtom().GetAtomicNum() == 16 else bond.GetEndAtom()
        if any(neighbor.GetAtomicNum() == 16 for neighbor in sulfur.GetNeighbors()):
            return None
        return "thioether"
    if numbers == [6, 6] and bond.GetBondType() == Chem.BondType.DOUBLE:
        return "olefin_metathesis"
    return None


def _partition_after_cuts(
    molecule: Chem.Mol,
    chain: Sequence[BackboneUnit],
    cut_bonds: set[frozenset[int]],
) -> tuple[list[set[int]], dict[int, int]] | None:
    """A valid decomposition has exactly one backbone Cα in every component."""
    components: list[set[int]] = []
    owner: dict[int, int] = {}
    for atom in molecule.GetAtoms():
        if atom.GetIdx() in owner:
            continue
        component: set[int] = set()
        queue = [atom.GetIdx()]
        while queue:
            current = queue.pop()
            if current in component:
                continue
            component.add(current)
            for neighbor in molecule.GetAtomWithIdx(current).GetNeighbors():
                if frozenset((current, neighbor.GetIdx())) not in cut_bonds:
                    queue.append(neighbor.GetIdx())
        anchors = [position for position, unit in enumerate(chain) if unit.ca in component]
        if len(anchors) != 1:
            return None
        position = anchors[0]
        components.append(component)
        owner.update({index: position for index in component})
    if len(components) != len(chain):
        return None
    return components, owner


def _crosslink_solution(
    molecule: Chem.Mol,
    chain: Sequence[BackboneUnit],
    backbone_kind: str,
) -> tuple[set[frozenset[int]], list[Chem.Bond], dict[int, int]] | None:
    peptide_cuts = {
        frozenset((chain[index].carbonyl, chain[index + 1].n))
        for index in range(len(chain) - 1)
    }
    if backbone_kind == "head_to_tail":
        peptide_cuts.add(frozenset((chain[-1].carbonyl, chain[0].n)))
    candidate_bonds = [
        bond for bond in molecule.GetBonds()
        if bond.IsInRing()
        and frozenset((bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())) not in peptide_cuts
        and _link_kind(bond) is not None
    ]
    # Bound the exact search. A molecule needing a larger search remains QC-visible
    # instead of silently receiving an arbitrary topology.
    if len(candidate_bonds) > 14:
        return None
    for count in range(len(candidate_bonds) + 1):
        valid: list[tuple[int, set[frozenset[int]], list[Chem.Bond], dict[int, int]]] = []
        for chosen in combinations(candidate_bonds, count):
            cuts = peptide_cuts | {
                frozenset((bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())) for bond in chosen
            }
            partition = _partition_after_cuts(molecule, chain, cuts)
            if partition is None:
                continue
            _, owner = partition
            if any(owner[bond.GetBeginAtomIdx()] == owner[bond.GetEndAtomIdx()] for bond in chosen):
                continue
            score = 0
            for bond in chosen:
                if _link_kind(bond) == "thioether":
                    sulfur = next(
                        atom.GetIdx() for atom in (bond.GetBeginAtom(), bond.GetEndAtom())
                        if atom.GetAtomicNum() == 16
                    )
                    path = _walk_path_without_cuts(molecule, chain[owner[sulfur]].ca, sulfur, cuts)
                    score += len(path) - 1 if path else 10000
            valid.append((score, cuts, list(chosen), owner))
        if valid:
            valid.sort(key=lambda item: item[0])
            if len(valid) > 1 and valid[0][0] == valid[1][0]:
                return None
            _, cuts, selected, owner = valid[0]
            return cuts, selected, owner
    return None


def _walk_path_without_cuts(
    molecule: Chem.Mol, start: int, finish: int, cuts: set[frozenset[int]]
) -> list[int]:
    queue = deque([(start, [start])])
    visited = {start}
    while queue:
        current, path = queue.popleft()
        if current == finish:
            return path
        for neighbor in molecule.GetAtomWithIdx(current).GetNeighbors():
            index = neighbor.GetIdx()
            if index not in visited and frozenset((current, index)) not in cuts:
                visited.add(index)
                queue.append((index, [*path, index]))
    return []


def _walk_component(
    molecule: Chem.Mol,
    seeds: Iterable[int],
    blocked: set[int],
    cut_bonds: set[frozenset[int]] | None = None,
) -> set[int]:
    stack = list(seeds)
    visited: set[int] = set()
    while stack:
        index = stack.pop()
        if index in blocked or index in visited:
            continue
        visited.add(index)
        atom = molecule.GetAtomWithIdx(index)
        stack.extend(
            neighbor.GetIdx() for neighbor in atom.GetNeighbors()
            if frozenset((index, neighbor.GetIdx())) not in (cut_bonds or set())
        )
    return visited


def _fragment_smiles(molecule: Chem.Mol, atom_indices: Iterable[int]) -> str | None:
    atoms = sorted(set(atom_indices))
    if not atoms:
        return None
    return Chem.MolFragmentToSmiles(
        molecule,
        atomsToUse=atoms,
        canonical=True,
        isomericSmiles=True,
        allBondsExplicit=False,
        allHsExplicit=False,
    )


def _external_neighbors(molecule: Chem.Mol, anchor: int, excluded: set[int]) -> list[int]:
    atom = molecule.GetAtomWithIdx(anchor)
    return sorted(neighbor.GetIdx() for neighbor in atom.GetNeighbors() if neighbor.GetIdx() not in excluded)


def _is_carbonyl_carbon(molecule: Chem.Mol, index: int) -> bool:
    atom = molecule.GetAtomWithIdx(index)
    return atom.GetAtomicNum() == 6 and bool(_carbonyl_oxygens(molecule, index))


def _cip_label(molecule: Chem.Mol, unit: BackboneUnit, side_atoms: set[int]) -> str:
    ca = molecule.GetAtomWithIdx(unit.ca)
    heavy_side_neighbors = [
        atom.GetIdx()
        for atom in ca.GetNeighbors()
        if atom.GetIdx() not in {unit.n, unit.carbonyl} and atom.GetAtomicNum() > 1
    ]
    if not heavy_side_neighbors:
        return "achiral"
    if len(heavy_side_neighbors) != 1:
        return "unknown"
    if not ca.HasProp("_CIPCode"):
        return "unknown"
    cip = ca.GetProp("_CIPCode")
    # L-cysteine is R because sulfur changes CIP priority; other supported
    # proteinogenic alpha-amino acids are L=S. D/L is never inferred from a
    # drawing orientation.
    first = molecule.GetAtomWithIdx(heavy_side_neighbors[0])
    cysteine_like = any(neighbor.GetAtomicNum() in (16, 34) for neighbor in first.GetNeighbors())
    if cysteine_like:
        return "L" if cip == "R" else "D" if cip == "S" else "unknown"
    return "L" if cip == "S" else "D" if cip == "R" else "unknown"


def _decompose_molecule(molecule: Chem.Mol, canonical: str, structure_hash: str) -> Decomposition:
    fragments = Chem.GetMolFrags(molecule)
    if len(fragments) != 1:
        return Decomposition(
            status="unsupported_multiple_fragments",
            canonical_smiles=canonical,
            structure_hash=structure_hash,
            length=None,
            n_cap=None,
            c_cap=None,
            residues=[],
            error="Multiple fragments are retained but not guessed into a main structure.",
            molecular_weight=Descriptors.MolWt(molecule),
        )

    units = _candidate_units(molecule)
    if not units:
        return Decomposition(
            status="excluded_non_peptide",
            canonical_smiles=canonical,
            structure_hash=structure_hash,
            length=0,
            n_cap=None,
            c_cap=None,
            residues=[],
            molecular_weight=Descriptors.MolWt(molecule),
        )

    outgoing, incoming = _unit_graph(molecule, units)
    paths = _maximal_paths(outgoing, incoming)
    cycles = _backbone_cycles(outgoing, units)
    longest_path = len(paths[0]) if paths else 0
    longest_cycle = len(cycles[0]) if cycles else 0
    backbone_kind = "head_to_tail" if longest_cycle > longest_path else "linear"
    candidates = cycles if backbone_kind == "head_to_tail" else paths
    if not candidates:
        return Decomposition(
            status="excluded_non_linear",
            canonical_smiles=canonical,
            structure_hash=structure_hash,
            length=None,
            n_cap=None,
            c_cap=None,
            residues=[],
            error="Candidate backbone is cyclic or has no unique N terminus.",
            molecular_weight=Descriptors.MolWt(molecule),
        )
    maximum = len(candidates[0])
    longest = [path for path in candidates if len(path) == maximum]
    atom_paths = {
        tuple(sorted((units[index].n, units[index].ca, units[index].carbonyl) for index in path))
        if backbone_kind == "head_to_tail"
        else tuple((units[index].n, units[index].ca, units[index].carbonyl) for index in path)
        for path in longest
    }
    if len(atom_paths) != 1:
        return Decomposition(
            status="unsupported_backbone",
            canonical_smiles=canonical,
            structure_hash=structure_hash,
            length=maximum,
            n_cap=None,
            c_cap=None,
            residues=[],
            error="Multiple equally long backbone candidates; no arbitrary choice was made.",
            molecular_weight=Descriptors.MolWt(molecule),
        )

    path = longest[0]
    if backbone_kind == "head_to_tail":
        start = min(range(len(path)), key=lambda index: units[path[index]].n)
        path = path[start:] + path[:start]
    chain = [units[index] for index in path]
    selected = set(path)
    if any(
        len([item for item in outgoing[index] if item in selected]) > 1
        or len([item for item in incoming[index] if item in selected]) > 1
        for index in path
    ):
        return Decomposition(
            status="excluded_non_linear",
            canonical_smiles=canonical,
            structure_hash=structure_hash,
            length=maximum,
            n_cap=None,
            c_cap=None,
            residues=[],
            error="Branched peptide backbone candidate.",
            molecular_weight=Descriptors.MolWt(molecule),
        )

    closure = molecule.GetBondBetweenAtoms(chain[-1].carbonyl, chain[0].n)
    if closure is not None and backbone_kind != "head_to_tail":
        return Decomposition(
            status="excluded_non_linear",
            canonical_smiles=canonical,
            structure_hash=structure_hash,
            length=maximum,
            n_cap=None,
            c_cap=None,
            residues=[],
            error="Cyclic backbone is outside the prototype scope.",
            molecular_weight=Descriptors.MolWt(molecule),
        )

    solution = _crosslink_solution(molecule, chain, backbone_kind)
    if solution is None:
        return Decomposition(
            status="unsupported_backbone",
            canonical_smiles=canonical,
            structure_hash=structure_hash,
            length=maximum,
            n_cap=None,
            c_cap=None,
            residues=[],
            error="Cyclization bonds cannot be cut into a unique set of residue fragments.",
            molecular_weight=Descriptors.MolWt(molecule),
        )
    cut_bonds, crosslinks, owner = solution
    # The virtual cuts must be exactly the inter-fragment bonds. Keeping the
    # original atom indices and bond records makes reconstruction lossless.
    original_edges = {
        frozenset((bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()))
        for bond in molecule.GetBonds()
    }
    inter_fragment_edges = {
        frozenset((bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()))
        for bond in molecule.GetBonds()
        if owner[bond.GetBeginAtomIdx()] != owner[bond.GetEndAtomIdx()]
    }
    if set(owner) != set(range(molecule.GetNumAtoms())) or inter_fragment_edges != cut_bonds:
        return Decomposition(
            status="unsupported_backbone", canonical_smiles=canonical,
            structure_hash=structure_hash, length=maximum, n_cap=None,
            c_cap=None, residues=[],
            error="Virtual cuts do not reconstruct the original atom and bond graph.",
            molecular_weight=Descriptors.MolWt(molecule),
        )
    all_anchors = {
        value
        for unit in chain
        for value in (unit.n, unit.ca, unit.carbonyl, *unit.oxygens)
    }
    n_cap_atoms: set[int] = set()
    c_cap_atoms: set[int] = set()
    if backbone_kind == "head_to_tail":
        n_cap = None
        c_cap = None
    else:
        n_neighbors = _external_neighbors(molecule, chain[0].n, {chain[0].ca})
        n_cap_seeds = [index for index in n_neighbors if _is_carbonyl_carbon(molecule, index)]
        n_cap_atoms = _walk_component(molecule, n_cap_seeds, all_anchors, cut_bonds)
        n_cap = _fragment_smiles(molecule, n_cap_atoms) or "H"
        last = chain[-1]
        c_cap_seeds = _external_neighbors(molecule, last.carbonyl, {last.ca, *last.oxygens})
        c_cap_atoms = _walk_component(molecule, c_cap_seeds, all_anchors, cut_bonds)
        c_cap = _fragment_smiles(molecule, c_cap_atoms) or "unspecified"

    ownership: list[set[int]] = []
    for unit in chain:
        stack = [unit.ca]
        atoms: set[int] = set()
        while stack:
            current = stack.pop()
            if current in atoms:
                continue
            atoms.add(current)
            current_atom = molecule.GetAtomWithIdx(current)
            for neighbor in current_atom.GetNeighbors():
                edge = frozenset((current, neighbor.GetIdx()))
                if edge not in cut_bonds:
                    stack.append(neighbor.GetIdx())
        ownership.append(atoms - n_cap_atoms - c_cap_atoms)

    # Cutting peptide bonds must partition the connected molecule into exactly
    # one component per residue (with terminal caps removed).
    for left in range(len(ownership)):
        for right in range(left + 1, len(ownership)):
            if ownership[left] & ownership[right]:
                return Decomposition(
                    status="unsupported_backbone",
                    canonical_smiles=canonical,
                    structure_hash=structure_hash,
                    length=maximum,
                    n_cap=n_cap,
                    c_cap=c_cap,
                    residues=[],
                    error="Backbone cuts did not produce disjoint residue atom sets.",
                    molecular_weight=Descriptors.MolWt(molecule),
                )

    residues: list[Residue] = []
    for position, (unit, atom_set) in enumerate(zip(chain, ownership, strict=True), start=1):
        side_seeds = _external_neighbors(molecule, unit.ca, {unit.n, unit.carbonyl})
        side_atoms = _walk_component(molecule, side_seeds, all_anchors, cut_bonds) & atom_set
        previous_carbonyl = chain[position - 2].carbonyl if position > 1 else None
        n_exclusions = {unit.ca}
        if previous_carbonyl is not None:
            n_exclusions.add(previous_carbonyl)
        r3_seeds = [
            index
            for index in _external_neighbors(molecule, unit.n, n_exclusions)
            if index not in n_cap_atoms
        ]
        r3_atoms = _walk_component(molecule, r3_seeds, all_anchors, cut_bonds) & atom_set
        bridge_atoms = side_atoms & r3_atoms
        bridge = _fragment_smiles(molecule, bridge_atoms)
        stereo = _cip_label(molecule, unit, side_atoms)
        side = _fragment_smiles(molecule, side_atoms)
        r3 = _fragment_smiles(molecule, r3_atoms)
        if bridge:
            residue_type = "Pro-like bridge"
            r1 = None
            r2 = None
            r3 = None
        elif stereo == "achiral" and r3:
            residue_type = "N-substituted Gly / peptoid"
            r1 = None
            r2 = None
        elif stereo == "achiral":
            residue_type = "Gly"
            r1 = None
            r2 = None
        elif len(side_seeds) > 1:
            residue_type = "alpha,alpha-disubstituted"
            r1 = side
            r2 = None
        elif r3:
            residue_type = "N-substituted amino acid"
            r1 = side if stereo != "D" else None
            r2 = side if stereo == "D" else None
        else:
            residue_type = "amino acid"
            r1 = side if stereo != "D" else None
            r2 = side if stereo == "D" else None
        variant_smiles = _fragment_smiles(molecule, atom_set) or ""
        variant_payload = json.dumps(
            {
                "fragment": variant_smiles,
                "stereo": stereo,
                "r1": r1,
                "r2": r2,
                "r3": r3,
                "bridge": bridge,
            },
            sort_keys=True,
            separators=(",", ":"),
        )
        residues.append(
            Residue(
                position=position,
                residue_type=residue_type,
                r1=r1,
                r2=r2,
                r3=r3,
                bridge=bridge,
                stereo_status=stereo,
                atom_indices=sorted(atom_set),
                anchors={
                    "n": unit.n,
                    "ca": unit.ca,
                    "carbonyl": unit.carbonyl,
                    "oxygens": list(unit.oxygens),
                },
                variant_hash=sha256_text(variant_payload),
            )
        )

    def endpoint(atom_index: int) -> dict[str, int | str]:
        if atom_index in n_cap_atoms:
            return {"kind": "n_cap", "position": 0, "atom_index": atom_index}
        if atom_index in c_cap_atoms:
            return {"kind": "c_cap", "position": len(chain) + 1, "atom_index": atom_index}
        return {"kind": "residue", "position": owner[atom_index] + 1, "atom_index": atom_index}

    links: list[dict] = []
    if backbone_kind == "head_to_tail":
        bond = molecule.GetBondBetweenAtoms(chain[-1].carbonyl, chain[0].n)
        links.append({
            "kind": "head_to_tail",
            "first": endpoint(chain[-1].carbonyl),
            "second": endpoint(chain[0].n),
            "bond_indices": [chain[-1].carbonyl, chain[0].n],
            "bond_order": 1.0,
            "stereo": "STEREONONE",
            "linker_smiles": "peptide_amide",
        })
    for bond in crosslinks:
        first_index, second_index = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        first_owner, second_owner = owner[first_index], owner[second_index]
        first_path = _walk_path_without_cuts(molecule, chain[first_owner].ca, first_index, cut_bonds)
        second_path = _walk_path_without_cuts(molecule, chain[second_owner].ca, second_index, cut_bonds)
        linker_atoms = set(first_path + second_path)
        if _link_kind(bond) == "triazole":
            # The shortest path crosses only part of the aromatic ring. Keep
            # its complete chemistry in the topology signature.
            linker_atoms.update(
                index
                for match in molecule.GetSubstructMatches(_TRIAZOLE_QUERY)
                if first_index == match[0] or second_index == match[0]
                for index in match
            )
        links.append({
            "kind": _link_kind(bond),
            "first": endpoint(first_index),
            "second": endpoint(second_index),
            "bond_indices": [first_index, second_index],
            "bond_order": bond.GetBondTypeAsDouble(),
            "stereo": str(bond.GetStereo()),
            "linker_smiles": _fragment_smiles(molecule, linker_atoms) or "",
        })

    def rotated_links(offset: int) -> list[tuple]:
        normalized: list[tuple] = []
        for link in links:
            ends = []
            for name in ("first", "second"):
                end = link[name]
                position = int(end["position"])
                if end["kind"] == "residue" and backbone_kind == "head_to_tail":
                    position = ((position - offset - 1) % len(residues)) + 1
                ends.append((end["kind"], position))
            normalized.append((
                link["kind"], tuple(sorted(ends)), link["bond_order"],
                link["stereo"], link["linker_smiles"],
            ))
        return sorted(normalized)

    canonical_start = 1
    if backbone_kind == "head_to_tail":
        offsets = range(len(residues))
        canonical_offset = min(
            offsets,
            key=lambda offset: (
                tuple(residues[(offset + index) % len(residues)].variant_hash for index in range(len(residues))),
                rotated_links(offset),
            ),
        )
        canonical_start = canonical_offset + 1
        residues = residues[canonical_offset:] + residues[:canonical_offset]
        for position, residue in enumerate(residues, start=1):
            residue.position = position
        for link in links:
            for name in ("first", "second"):
                end = link[name]
                if end["kind"] == "residue":
                    end["position"] = ((int(end["position"]) - canonical_offset - 1) % len(residues)) + 1
    signature = sha256_text(json.dumps({
        "backbone_kind": backbone_kind,
        "length": len(residues),
        "links": min((rotated_links(offset) for offset in range(len(residues))), default=[])
        if backbone_kind == "head_to_tail" else rotated_links(0),
    }, sort_keys=True, separators=(",", ":")))
    status = "accepted" if len(chain) >= 5 else "excluded_short_peptide"
    return Decomposition(
        status=status,
        canonical_smiles=canonical,
        structure_hash=structure_hash,
        length=len(chain),
        n_cap=n_cap,
        c_cap=c_cap,
        residues=residues,
        molecular_weight=Descriptors.MolWt(molecule),
        backbone_kind=backbone_kind,
        cyclization_links=links,
        topology_signature=signature,
        canonical_start=canonical_start,
    )


def decompose_smiles(smiles: str) -> Decomposition:
    try:
        molecule, canonical, structure_hash = canonicalize_smiles(smiles)
    except ValueError as exc:
        return Decomposition(
            status="invalid_structure",
            canonical_smiles=None,
            structure_hash=None,
            length=None,
            n_cap=None,
            c_cap=None,
            residues=[],
            error=str(exc),
        )
    # Persisted atom indices must address the canonical SMILES served by the API,
    # regardless of the input SMILES traversal order.
    molecule = Chem.MolFromSmiles(canonical)
    if molecule is None:
        raise ValueError("Canonical SMILES could not be parsed")
    Chem.AssignStereochemistry(molecule, cleanIt=True, force=True)
    try:
        return _decompose_molecule(molecule, canonical, structure_hash)
    except Exception as exc:  # a single unusual molecule must not abort an import
        return Decomposition(
            status="unsupported_backbone",
            canonical_smiles=canonical,
            structure_hash=structure_hash,
            length=None,
            n_cap=None,
            c_cap=None,
            residues=[],
            error=f"Conservative parser failure: {type(exc).__name__}: {exc}",
            molecular_weight=Descriptors.MolWt(molecule),
        )


def formula(smiles: str) -> str:
    molecule = Chem.MolFromSmiles(smiles)
    return rdMolDescriptors.CalcMolFormula(molecule) if molecule is not None else ""


def structure_svg(
    smiles: str,
    highlight_atoms: Sequence[int] = (),
    width: int = 520,
    height: int = 260,
) -> str:
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        raise ValueError("Invalid molecule")
    rdDepictor.Compute2DCoords(molecule)
    drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
    options = drawer.drawOptions()
    options.addStereoAnnotation = True
    drawer.DrawMolecule(molecule, highlightAtoms=list(highlight_atoms))
    drawer.FinishDrawing()
    return drawer.GetDrawingText()


def _linear_scaffold_molecule(smiles: str, residues: Sequence[Residue]) -> Chem.Mol:
    """Constrain every backbone bond to an N-to-C zigzag while laying out side chains."""
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None or not residues:
        raise ValueError("Accepted peptide structure and residues are required")
    anchors = {}
    for index, residue in enumerate(residues):
        for offset, anchor in enumerate(("n", "ca", "carbonyl")):
            backbone_index = 3 * index + offset
            anchors[int(residue.anchors[anchor])] = Point2D(
                backbone_index * 1.3,
                0.75 if backbone_index % 2 else 0,
            )
    rdDepictor.Compute2DCoords(molecule, canonOrient=False, coordMap=anchors)
    return molecule


def _cyclic_scaffold_molecule(smiles: str, residues: Sequence[Residue]) -> Chem.Mol:
    """Let CoordGen place the complete macrocycle without distorting bonds."""
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None or not residues:
        raise ValueError("Accepted cyclic peptide structure and residues are required")
    rdCoordGen.AddCoords(molecule)
    return molecule


def parent_scaffold_svg(
    smiles: str,
    residues: Sequence[Residue],
    width: int = 640,
    height: int = 220,
    changed_position: int | None = None,
    site_kind: str = "residue",
    cyclization_links: Sequence[dict] | None = None,
) -> str:
    """Draw a peptide scaffold with the changed block in yellow."""
    if cyclization_links:
        molecule = _cyclic_scaffold_molecule(smiles, residues)
        for residue in residues:
            ca = int(residue.anchors["ca"])
            if 0 <= ca < molecule.GetNumAtoms():
                molecule.GetAtomWithIdx(ca).SetProp("atomNote", str(residue.position))
        highlighted = (
            set(residues[changed_position - 1].atom_indices)
            if site_kind == "residue" and changed_position is not None and 1 <= changed_position <= len(residues)
            else set()
        )
        yellow = (1.0, 0.88, 0.2)
        palette = {
            "head_to_tail": (0.27, 0.48, 0.83),
            "disulfide": (0.55, 0.36, 0.72),
            "thioether": (0.14, 0.62, 0.53),
            "olefin_metathesis": (0.90, 0.42, 0.28),
            "triazole": (0.72, 0.40, 0.70),
        }
        bond_colors: dict[int, tuple[float, float, float]] = {}
        for link in cyclization_links:
            atoms = link.get("bond_indices", [])
            if len(atoms) != 2:
                continue
            bond = molecule.GetBondBetweenAtoms(int(atoms[0]), int(atoms[1]))
            if bond is not None:
                bond_colors[bond.GetIdx()] = palette.get(link["kind"], (0.5, 0.5, 0.5))
        for bond in molecule.GetBonds():
            if bond.GetBeginAtomIdx() in highlighted and bond.GetEndAtomIdx() in highlighted:
                bond_colors[bond.GetIdx()] = yellow
        Chem.WedgeMolBonds(molecule, molecule.GetConformer())
        drawing_side = min(900, max(420, len(residues) * 40))
        drawing_width = max(width, drawing_side)
        drawing_height = max(height, drawing_side)
        drawer = rdMolDraw2D.MolDraw2DSVG(drawing_width, drawing_height)
        drawer.drawOptions().addStereoAnnotation = True
        drawer.DrawMolecule(
            molecule,
            highlightAtoms=sorted(highlighted),
            highlightBonds=sorted(bond_colors),
            highlightAtomColors={index: yellow for index in highlighted},
            highlightBondColors=bond_colors,
        )
        drawer.FinishDrawing()
        return drawer.GetDrawingText()
    molecule = _linear_scaffold_molecule(smiles, residues)
    conformer = molecule.GetConformer()
    highlighted: set[int] = set()
    if site_kind == "residue" and changed_position is not None and 1 <= changed_position <= len(residues):
        highlighted.update(residues[changed_position - 1].atom_indices)
    elif site_kind in {"n_cap", "c_cap"}:
        core = {atom for residue in residues for atom in residue.atom_indices}
        terminal = int(residues[0].anchors["n"] if site_kind == "n_cap" else residues[-1].anchors["carbonyl"])
        pending = [neighbor.GetIdx() for neighbor in molecule.GetAtomWithIdx(terminal).GetNeighbors() if neighbor.GetIdx() not in core]
        while pending:
            atom = pending.pop()
            if atom in highlighted:
                continue
            highlighted.add(atom)
            pending.extend(neighbor.GetIdx() for neighbor in molecule.GetAtomWithIdx(atom).GetNeighbors() if neighbor.GetIdx() not in core and neighbor.GetIdx() not in highlighted)
        highlighted.add(terminal)
    Chem.WedgeMolBonds(molecule, conformer)
    drawer = rdMolDraw2D.MolDraw2DSVG(max(width, len(residues) * 72), height)
    drawer.drawOptions().addStereoAnnotation = True
    bonds = [bond.GetIdx() for bond in molecule.GetBonds() if bond.GetBeginAtomIdx() in highlighted and bond.GetEndAtomIdx() in highlighted]
    yellow = (1.0, 0.88, 0.2)
    drawer.DrawMolecule(
        molecule,
        highlightAtoms=sorted(highlighted),
        highlightBonds=bonds,
        highlightAtomColors={atom: yellow for atom in highlighted},
        highlightBondColors={bond: yellow for bond in bonds},
    )
    drawer.FinishDrawing()
    return drawer.GetDrawingText()


def residue_building_block_svg(
    smiles: str, residue: Residue, width: int = 300, height: int = 200,
    attachment_atoms: Sequence[int] = (),
) -> str:
    """Draw only the residue atom set, after cutting its peptide-bond connections."""
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        raise ValueError("Invalid peptide structure")
    fragment = _fragment_smiles(molecule, residue.atom_indices)
    if not fragment:
        raise ValueError("Empty residue fragment")
    if not attachment_atoms:
        return structure_svg(fragment, width=width, height=height)
    kept = sorted(set(residue.atom_indices))
    index_map = {old: new for new, old in enumerate(kept)}
    drawing = Chem.RWMol(molecule)
    for index in reversed(range(molecule.GetNumAtoms())):
        if index not in index_map:
            drawing.RemoveAtom(index)
    fragment_molecule = drawing.GetMol()
    for old_index in attachment_atoms:
        if old_index in index_map:
            fragment_molecule.GetAtomWithIdx(index_map[old_index]).SetProp("atomNote", "bridge *")
    rdDepictor.Compute2DCoords(fragment_molecule)
    drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
    drawer.drawOptions().addStereoAnnotation = True
    drawer.DrawMolecule(fragment_molecule)
    drawer.FinishDrawing()
    return drawer.GetDrawingText()


def aligned_pair_svgs(
    first_smiles: str,
    second_smiles: str,
    first_residues: Sequence[Residue],
    second_residues: Sequence[Residue],
    changed_position: int,
) -> tuple[str, str, str | None]:
    first = Chem.MolFromSmiles(first_smiles)
    second = Chem.MolFromSmiles(second_smiles)
    if first is None or second is None:
        raise ValueError("Invalid pair structure")
    rdDepictor.Compute2DCoords(first)
    conformer = first.GetConformer()
    first_n = int(first_residues[0].anchors["n"])
    last_c = int(first_residues[-1].anchors["carbonyl"])
    if conformer.GetAtomPosition(first_n).x > conformer.GetAtomPosition(last_c).x:
        for index in range(first.GetNumAtoms()):
            point = conformer.GetAtomPosition(index)
            point.x = -point.x
            point.y = -point.y
            conformer.SetAtomPosition(index, point)
    atom_map: list[tuple[int, int]] = []
    for left, right in zip(first_residues, second_residues, strict=True):
        for key in ("n", "ca", "carbonyl"):
            atom_map.append((int(right.anchors[key]), int(left.anchors[key])))
    warning: str | None = None
    try:
        rdDepictor.GenerateDepictionMatching2DStructure(second, first, atomMap=atom_map)
    except Exception as exc:
        rdDepictor.Compute2DCoords(second)
        warning = f"Backbone alignment unavailable; full structures shown ({type(exc).__name__})."
    Chem.WedgeMolBonds(first, first.GetConformer())
    Chem.WedgeMolBonds(second, second.GetConformer())

    def draw(molecule: Chem.Mol, residue: Residue) -> str:
        drawer = rdMolDraw2D.MolDraw2DSVG(560, 280)
        drawer.drawOptions().addStereoAnnotation = True
        drawer.DrawMolecule(molecule, highlightAtoms=residue.atom_indices)
        drawer.FinishDrawing()
        return drawer.GetDrawingText()

    return (
        draw(first, first_residues[changed_position - 1]),
        draw(second, second_residues[changed_position - 1]),
        warning,
    )

