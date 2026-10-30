from __future__ import annotations

import argparse
import hashlib
from pathlib import Path

from app.db import connect, initialize, utcnow
from app.project_names import unique_project_name
from app.service import stage_import
from app.worker import run_job


def main() -> None:
    parser = argparse.ArgumentParser(description="Load the bundled MDM2 and GLP-1 prototypes into separate local projects.")
    parser.add_argument("--file", type=Path, action="append", help="CSV to seed (repeatable; defaults to the dated MDM2 and GLP-1 sources in data/source)")
    parser.add_argument("--project-name", help="Filename or label used after the yyyymmdd_ prefix (only with one --file)")
    arguments = parser.parse_args()
    files = arguments.file or [Path("data/source/20261030_MDM2.csv"), Path("data/source/20261030_GLP-1.csv")]
    if arguments.project_name and len(files) != 1:
        parser.error("--project-name requires exactly one --file")
    for path in files:
        if not path.is_file():
            parser.error(f"input file not found: {path}")
        if path.suffix.lower() != ".csv":
            parser.error(f"CSV file required: {path}")
    initialize()
    with connect() as db:
        for path in files:
            content = path.read_bytes()
            digest = hashlib.sha256(content).hexdigest()
            project = db.execute(
                """SELECT p.project_id FROM projects p JOIN imports i ON i.project_id=p.project_id
                   WHERE i.source_name=? ORDER BY p.project_id, i.import_id LIMIT 1""",
                (path.name,),
            ).fetchone()
            if project is None:
                now = utcnow()
                name = unique_project_name(db, arguments.project_name or path.name, now)
                project_id = db.execute(
                    "INSERT INTO projects(name, settings_json, created_at, updated_at) VALUES (?, '{}', ?, ?)",
                    (name, now, now),
                ).lastrowid
            else:
                project_id = project["project_id"]
            already_complete = db.execute(
                """SELECT 1 FROM imports WHERE project_id=? AND file_sha256=? AND import_status='complete'
                   LIMIT 1""",
                (project_id, digest),
            ).fetchone()
            if already_complete:
                print(f"Project {project_id} ({path.name}): identical source already analyzed; skipped.")
                continue
            already_pending = db.execute(
                """SELECT j.status FROM imports i JOIN jobs j ON j.import_id=i.import_id
                   WHERE i.project_id=? AND i.file_sha256=? AND j.status IN ('queued', 'running')
                   ORDER BY i.import_id DESC LIMIT 1""",
                (project_id, digest),
            ).fetchone()
            if already_pending:
                print(f"Project {project_id} ({path.name}): identical source is already {already_pending['status']}; skipped.")
                continue
            staged = stage_import(db, project_id, path.name, "csv", content)
            if staged["duplicate"]:
                print(f"Project {project_id} ({path.name}): identical source already queued; skipped.")
                continue
            job_id = int(staged["job_id"])
            claimed = db.execute(
                """UPDATE jobs SET status='running', progress=0.01, updated_at=?
                   WHERE job_id=? AND status='queued'""",
                (utcnow(), job_id),
            ).rowcount
            if not claimed:
                print(f"Project {project_id} ({path.name}): queued job {job_id} was claimed by the worker.")
                continue
            run_job({"job_id": job_id, "job_type": "import_and_analyze", "import_id": staged["import_id"]})
            completed = db.execute("SELECT status, metrics_json, error FROM jobs WHERE job_id=?", (job_id,)).fetchone()
            if completed["status"] != "complete":
                raise RuntimeError(f"Project {project_id} ({path.name}) analysis failed: {completed['error']}")
            print(f"Project {project_id} ({path.name}): {completed['metrics_json']}")


if __name__ == "__main__":
    main()

