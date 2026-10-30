"""Stage local cyclic literature CSVs as one ring-Peptide snapshot."""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
from pathlib import Path

from app.db import connect, initialize, utcnow
from app.importers import read_csv
from app.service import stage_import


SOURCE_STEMS = ("Grb2-SH2_S", "MEGA-A4_S", "Oxytocin_S", "RGD_H2T", "UbE2")


def combined_csv(source_dir: Path, date: str) -> tuple[bytes, int]:
    source_paths = [source_dir / f"{date}_{stem}.csv" for stem in SOURCE_STEMS]
    rap1 = next((path for path in (source_dir / f"{date}_RAP1.csv", source_dir / "RAP1.csv") if path.is_file()), None)
    if rap1 is not None:
        source_paths.append(rap1)
    duplicate = source_dir / f"{date}_oxitosin_SS.csv"
    for path in source_paths:
        if not path.is_file():
            raise FileNotFoundError(path)
    if duplicate.is_file() and duplicate.read_bytes() != source_paths[2].read_bytes():
        raise ValueError(f"Expected duplicate has different content: {duplicate}")

    fields: list[str] | None = None
    rows: list[dict[str, str]] = []
    for path in source_paths:
        records = list(read_csv(path))
        if not records:
            raise ValueError(f"Source CSV has no records: {path}")
        if fields is None:
            fields = list(records[0].properties)
        for record in records:
            if list(record.properties) != fields:
                raise ValueError(f"CSV columns differ: {path}")
            rows.append({**record.properties, "Source File": path.name})

    output = io.StringIO(newline="")
    writer = csv.DictWriter(output, fieldnames=[*fields, "Source File"], delimiter=";", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return output.getvalue().encode("utf-8"), len(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-dir", type=Path, default=Path("data/source"))
    parser.add_argument("--date", default="20261030", help="yyyymmdd prefix on source filenames")
    args = parser.parse_args()
    if len(args.date) != 8 or not args.date.isdigit():
        parser.error("--date must be yyyymmdd")
    content, count = combined_csv(args.source_dir, args.date)
    digest = hashlib.sha256(content).hexdigest()
    initialize()
    with connect() as db:
        project = db.execute("SELECT project_id FROM projects WHERE name='ring-Peptide'").fetchone()
        if project is None:
            now = utcnow()
            project_id = db.execute(
                "INSERT INTO projects(name, created_at, updated_at) VALUES ('ring-Peptide', ?, ?)",
                (now, now),
            ).lastrowid
        else:
            project_id = project["project_id"]
        existing = db.execute(
            """SELECT import_status FROM imports WHERE project_id=? AND file_sha256=?
               AND import_status IN ('queued', 'running', 'complete')
               ORDER BY import_id DESC LIMIT 1""",
            (project_id, digest),
        ).fetchone()
        if existing:
            print(f"ring-Peptide: identical {count}-row source is already {existing['import_status']}; skipped.")
            return
        staged = stage_import(db, project_id, f"{args.date}_ring-Peptide.csv", "csv", content)
        print(f"ring-Peptide: staged {count} rows as import {staged['import_id']} (job {staged['job_id']}).")


if __name__ == "__main__":
    main()
