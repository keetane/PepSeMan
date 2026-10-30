"""Catch dependency license metadata changes before a publication update."""

from collections import Counter
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_javascript_lockfile_license_inventory_matches_review() -> None:
    lock = json.loads((ROOT / "frontend" / "package-lock.json").read_text())
    packages = {path: entry for path, entry in lock["packages"].items() if path}
    counts = Counter(entry.get("license", "MISSING") for entry in packages.values())
    assert counts == {
        "MIT": 282,
        "ISC": 32,
        "BSD-3-Clause": 21,
        "BSD-2-Clause": 6,
        "Apache-2.0": 6,
        "(MIT OR Apache-2.0)": 1,
        "CC-BY-4.0": 1,
        "BlueOak-1.0.0": 1,
        "Unlicense": 1,
        "Zlib": 1,
        "0BSD": 1,
        "MISSING": 2,
    }
    assert {
        path: entry["version"]
        for path, entry in packages.items()
        if "license" not in entry
    } == {
        "node_modules/glsl-resolve/node_modules/xtend": "2.2.0",
        "node_modules/stack-trace": "0.0.9",
    }
