from __future__ import annotations

import csv
import io
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator

from rdkit import Chem


@dataclass
class SourceRecord:
    row_number: int
    properties: dict[str, str]
    smiles: str
    source_molecule_id: str | None
    source_activity_id: str | None


ALIASES = {
    "smiles": ("Smiles", "SMILES", "smiles"),
    "molecule_id": ("Molecule ChEMBL ID", "Compound ChEMBL ID", "molecule_id", "ID", "id"),
    "activity_id": ("Activity ChEMBL ID", "activity_id"),
}


def _first(properties: dict[str, str], aliases: tuple[str, ...]) -> str | None:
    for name in aliases:
        value = properties.get(name)
        if value not in (None, "", "None"):
            return value
    return None


def read_csv(path: Path) -> Iterator[SourceRecord]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        sample = handle.read(8192)
        handle.seek(0)
        try:
            dialect = csv.Sniffer().sniff(sample, delimiters=";,\t,")
        except csv.Error:
            dialect = csv.excel
            dialect.delimiter = ";"
        reader = csv.DictReader(handle, dialect=dialect)
        if not reader.fieldnames:
            raise ValueError("CSV has no header")
        smiles_column = next((name for name in ALIASES["smiles"] if name in reader.fieldnames), None)
        if smiles_column is None:
            raise ValueError("CSV requires a Smiles/SMILES column")
        for row_number, row in enumerate(reader, start=2):
            properties = {str(key): value or "" for key, value in row.items() if key is not None}
            yield SourceRecord(
                row_number=row_number,
                properties=properties,
                smiles=properties[smiles_column],
                source_molecule_id=_first(properties, ALIASES["molecule_id"]),
                source_activity_id=_first(properties, ALIASES["activity_id"]),
            )


def read_smi(path: Path) -> Iterator[SourceRecord]:
    with path.open("r", encoding="utf-8") as handle:
        for row_number, line in enumerate(handle, start=1):
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            parts = stripped.split()
            smiles = parts[0]
            identifier = parts[1] if len(parts) > 1 else f"row-{row_number}"
            yield SourceRecord(
                row_number=row_number,
                properties={"Smiles": smiles, "ID": identifier},
                smiles=smiles,
                source_molecule_id=identifier,
                source_activity_id=None,
            )


def read_sdf(path: Path) -> Iterator[SourceRecord]:
    supplier = Chem.SDMolSupplier(str(path), removeHs=False, sanitize=True)
    for zero_index, molecule in enumerate(supplier):
        row_number = zero_index + 1
        if molecule is None:
            yield SourceRecord(
                row_number=row_number,
                properties={"_ParseError": "RDKit failed to parse SDF record"},
                smiles="",
                source_molecule_id=f"sdf-row-{row_number}",
                source_activity_id=None,
            )
            continue
        properties = {name: molecule.GetProp(name) for name in molecule.GetPropNames(includePrivate=False)}
        smiles = Chem.MolToSmiles(molecule, canonical=False, isomericSmiles=True)
        properties.setdefault("Smiles", smiles)
        yield SourceRecord(
            row_number=row_number,
            properties=properties,
            smiles=smiles,
            source_molecule_id=_first(properties, ALIASES["molecule_id"]) or f"sdf-row-{row_number}",
            source_activity_id=_first(properties, ALIASES["activity_id"]),
        )


def iter_records(path: Path, source_format: str) -> Iterator[SourceRecord]:
    if source_format == "csv":
        yield from read_csv(path)
    elif source_format == "smi":
        yield from read_smi(path)
    elif source_format == "sdf":
        yield from read_sdf(path)
    else:
        raise ValueError(f"Unsupported source format: {source_format}")


def csv_to_sdf(csv_path: Path, sdf_path: Path) -> int:
    writer = Chem.SDWriter(str(sdf_path))
    count = 0
    try:
        for record in read_csv(csv_path):
            molecule = Chem.MolFromSmiles(record.smiles)
            if molecule is None:
                continue
            from rdkit.Chem import rdDepictor

            rdDepictor.Compute2DCoords(molecule)
            for key, value in record.properties.items():
                if value:
                    molecule.SetProp(key, value)
            molecule.SetProp("_OriginalRowNumber", str(record.row_number))
            writer.write(molecule)
            count += 1
    finally:
        writer.close()
    return count


def semantic_roundtrip_report(csv_path: Path, sdf_path: Path) -> dict[str, int | list[str]]:
    csv_records = list(read_csv(csv_path))
    sdf_records = list(read_sdf(sdf_path))
    errors: list[str] = []
    if len(csv_records) != len(sdf_records):
        errors.append(f"record count differs: CSV={len(csv_records)} SDF={len(sdf_records)}")
    for index, (left, right) in enumerate(zip(csv_records, sdf_records, strict=False), start=1):
        left_mol = Chem.MolFromSmiles(left.smiles)
        right_mol = Chem.MolFromSmiles(right.smiles)
        if left_mol is None or right_mol is None:
            errors.append(f"record {index}: parse failure")
            continue
        left_canonical = Chem.MolToSmiles(left_mol, canonical=True, isomericSmiles=True)
        right_canonical = Chem.MolToSmiles(right_mol, canonical=True, isomericSmiles=True)
        if left_canonical != right_canonical:
            errors.append(f"record {index}: canonical isomeric SMILES differs")
        for key, value in left.properties.items():
            if value and right.properties.get(key) != value:
                errors.append(f"record {index}: property {key!r} differs")
                break
    return {"csv_records": len(csv_records), "sdf_records": len(sdf_records), "errors": errors}

