from __future__ import annotations

import hashlib
import json
import sqlite3
from typing import Any

from app.config import MMP_ENGINE_VERSION, PARSER_VERSION
from app.db import json_text, utcnow


def _hash(payload: Any) -> str:
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest()


def _sequence(db: sqlite3.Connection, structure_hash: str) -> dict[str, Any]:
    decomposition = db.execute(
        """SELECT n_cap, c_cap, backbone_kind, topology_signature, cyclization_links_json
           FROM decompositions
           WHERE structure_hash=? AND parser_version=? AND status='accepted'""",
        (structure_hash, PARSER_VERSION),
    ).fetchone()
    if decomposition is None:
        raise ValueError(f"Accepted decomposition not found: {structure_hash}")
    residues = db.execute(
        """SELECT variant_hash FROM residues
           WHERE structure_hash=? AND parser_version=? ORDER BY position""",
        (structure_hash, PARSER_VERSION),
    ).fetchall()
    return {
        "n_cap": decomposition["n_cap"],
        "residues": [row["variant_hash"] for row in residues],
        "c_cap": decomposition["c_cap"],
        "backbone_kind": decomposition["backbone_kind"],
        "topology_signature": decomposition["topology_signature"],
        "links": json.loads(decomposition["cyclization_links_json"]),
    }


def _link_signature(sequence: dict[str, Any], offset: int = 0) -> list[tuple]:
    length = len(sequence["residues"])
    result = []
    for link in sequence["links"]:
        if link["kind"] == "head_to_tail":
            continue
        ends = []
        for endpoint in (link["first"], link["second"]):
            position = endpoint["position"]
            if sequence["backbone_kind"] == "head_to_tail" and endpoint["kind"] == "residue":
                position = ((position - offset - 1) % length) + 1
            ends.append((endpoint["kind"], position))
        result.append((link["kind"], tuple(sorted(ends)), link["bond_order"], link["stereo"], link["linker_smiles"]))
    return sorted(result)


def _matching_sites(first: dict[str, Any], second: dict[str, Any]) -> list[tuple[str, int, int]]:
    if first["backbone_kind"] != second["backbone_kind"] or first["topology_signature"] != second["topology_signature"]:
        return []
    a, b = first["residues"], second["residues"]
    if len(a) != len(b):
        return []
    if first["backbone_kind"] == "linear":
        if _link_signature(first) != _link_signature(second):
            return []
        tokens_a = [first["n_cap"], *a, first["c_cap"]]
        tokens_b = [second["n_cap"], *b, second["c_cap"]]
        differences = [index for index, (left, right) in enumerate(zip(tokens_a, tokens_b, strict=True)) if left != right]
        if len(differences) != 1:
            return []
        index = differences[0]
        kind = "n_cap" if index == 0 else "c_cap" if index == len(tokens_a) - 1 else "residue"
        position = 0 if kind == "n_cap" else len(a) + 1 if kind == "c_cap" else index
        return [(kind, position, position)]
    matches: set[tuple[str, int, int]] = set()
    for shift in range(len(a)):
        if _link_signature(first) != _link_signature(second, shift):
            continue
        rotated = b[shift:] + b[:shift]
        differences = [index for index, (left, right) in enumerate(zip(a, rotated, strict=True)) if left != right]
        if len(differences) == 1:
            first_position = differences[0] + 1
            second_position = ((differences[0] + shift) % len(a)) + 1
            matches.add(("residue", first_position, second_position))
    return list(matches) if len(matches) == 1 else []


def context_rows(db: sqlite3.Connection, project_id: int, structure_hash: str) -> list[tuple]:
    sequence = _sequence(db, structure_hash)
    residues = sequence["residues"]
    cyclic = sequence["backbone_kind"] == "head_to_tail"
    tokens = [*(f"R:{item}" for item in residues)] if cyclic else [
        f"N:{sequence['n_cap']}", *(f"R:{item}" for item in residues), f"C:{sequence['c_cap']}"
    ]
    rows: list[tuple] = []
    for index, variant in enumerate(tokens):
        masked = [*tokens]
        masked[index] = "*"
        site_kind = "residue" if cyclic else "n_cap" if index == 0 else "c_cap" if index == len(tokens) - 1 else "residue"
        site_position = index + 1 if cyclic else 0 if index == 0 else len(residues) + 1 if index == len(tokens) - 1 else index
        variants = [
            (masked[offset:] + masked[:offset], _link_signature(sequence, offset))
            for offset in range(len(residues))
        ] if cyclic else [(masked, _link_signature(sequence))]
        rows.append(
            (
                project_id,
                structure_hash,
                site_kind,
                site_position,
                _hash((sequence["topology_signature"], min(variants))),
                _hash(variant),
                MMP_ENGINE_VERSION,
            )
        )
    return rows


def update_project_contexts(
    db: sqlite3.Connection,
    project_id: int,
    active_hashes: set[str],
) -> dict[str, int]:
    existing = {
        row["structure_hash"]
        for row in db.execute(
            "SELECT DISTINCT structure_hash FROM contexts WHERE project_id=? AND engine_version=?",
            (project_id, MMP_ENGINE_VERSION),
        )
    }
    added = active_hashes - existing
    removed = existing - active_hashes
    affected_contexts: set[str] = set()
    for structure_hash in removed:
        affected_contexts.update(
            row["context_hash"]
            for row in db.execute(
                "SELECT context_hash FROM contexts WHERE project_id=? AND structure_hash=? AND engine_version=?",
                (project_id, structure_hash, MMP_ENGINE_VERSION),
            )
        )
        db.execute(
            "DELETE FROM contexts WHERE project_id=? AND structure_hash=? AND engine_version=?",
            (project_id, structure_hash, MMP_ENGINE_VERSION),
        )
    for structure_hash in added:
        rows = context_rows(db, project_id, structure_hash)
        affected_contexts.update(row[4] for row in rows)
        db.executemany(
            """INSERT OR REPLACE INTO contexts(
                   project_id, structure_hash, site_kind, site_position,
                   context_hash, variant_hash, engine_version
               ) VALUES (?, ?, ?, ?, ?, ?, ?)""",
            rows,
        )

    created = 0
    sequence_cache: dict[str, dict[str, Any]] = {}

    def cached_sequence(structure_hash: str) -> dict[str, Any]:
        if structure_hash not in sequence_cache:
            sequence_cache[structure_hash] = _sequence(db, structure_hash)
        return sequence_cache[structure_hash]

    for context_hash in affected_contexts:
        db.execute(
            "DELETE FROM pairs WHERE project_id=? AND context_hash=? AND engine_version=?",
            (project_id, context_hash, MMP_ENGINE_VERSION),
        )
        bucket = db.execute(
            """SELECT structure_hash, site_kind, site_position, variant_hash
               FROM contexts
               WHERE project_id=? AND context_hash=? AND engine_version=?
               ORDER BY structure_hash""",
            (project_id, context_hash, MMP_ENGINE_VERSION),
        ).fetchall()
        for left_index, left in enumerate(bucket):
            for right in bucket[left_index + 1 :]:
                if left["variant_hash"] == right["variant_hash"]:
                    continue
                first_hash, second_hash = sorted((left["structure_hash"], right["structure_hash"]))
                first_seq, second_seq = cached_sequence(first_hash), cached_sequence(second_hash)
                sites = _matching_sites(first_seq, second_seq)
                if not sites:
                    continue
                site_kind, first_position, second_position = sites[0]
                pair_key = _hash((first_hash, second_hash, site_kind, first_position, second_position))
                inserted = db.execute(
                    """INSERT OR IGNORE INTO pairs(
                           project_id, pair_key, first_structure_hash, second_structure_hash,
                           site_kind, site_position, context_hash, engine_version, created_at,
                           second_site_position
                       ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)""",
                    (
                        project_id, pair_key, first_hash, second_hash, site_kind,
                        first_position, context_hash, MMP_ENGINE_VERSION, utcnow(),
                        second_position,
                    ),
                ).rowcount
                created += inserted
    return {
        "added_structures": len(added),
        "removed_structures": len(removed),
        "affected_context_buckets": len(affected_contexts),
        "pairs_created": created,
        "unchanged_structures": len(active_hashes & existing),
    }

