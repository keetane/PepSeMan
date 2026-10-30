"""Snapshot or verify the shared, non-natural-residue dictionary template."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from app.api import _dictionary_rows
from app.db import connect


TEMPLATE_PATH = Path(__file__).resolve().parents[1] / "templates" / "sequence_dictionary.json"
FIELDS = (
    "variant_hash",
    "residue_type",
    "r1",
    "r2",
    "r3",
    "bridge",
    "stereo_status",
    "default_abbreviation",
    "abbreviation",
    "custom_abbreviation",
)


def snapshot() -> dict[str, object]:
    with connect() as db:
        entries = [{key: row[key] for key in FIELDS} for row in _dictionary_rows(db)]
    return {"schema_version": 1, "sources": ["MDM2.csv", "GLP-1.csv"], "entries": entries}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true", help="Compare the current database with the tracked template")
    args = parser.parse_args()
    current = snapshot()
    if args.check:
        if json.loads(TEMPLATE_PATH.read_text(encoding="utf-8")) != current:
            parser.exit(1, "Dictionary differs from the tracked template.\n")
        print(f"Dictionary template matches {len(current['entries'])} building blocks.")
        return
    TEMPLATE_PATH.parent.mkdir(parents=True, exist_ok=True)
    TEMPLATE_PATH.write_text(json.dumps(current, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(f"Saved {len(current['entries'])} building blocks to {TEMPLATE_PATH}.")


if __name__ == "__main__":
    main()
