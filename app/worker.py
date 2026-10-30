from __future__ import annotations

import json
import signal
import time
from typing import Any

from app.db import connect, initialize, json_text, transaction, utcnow
from app.service import process_import


running = True


def _stop(_signum, _frame) -> None:
    global running
    running = False


def claim_job() -> dict[str, Any] | None:
    with connect() as db, transaction(db):
        row = db.execute(
            "SELECT * FROM jobs WHERE status='queued' ORDER BY created_at, job_id LIMIT 1"
        ).fetchone()
        if row is None:
            return None
        changed = db.execute(
            """UPDATE jobs SET status='running', progress=0.01, updated_at=?
               WHERE job_id=? AND status='queued'""",
            (utcnow(), row["job_id"]),
        ).rowcount
        return dict(row) if changed else None


def update_progress(job_id: int, value: float) -> None:
    with connect() as db:
        db.execute(
            "UPDATE jobs SET progress=?, updated_at=? WHERE job_id=? AND status='running'",
            (max(0.0, min(0.99, value)), utcnow(), job_id),
        )


def run_job(job: dict[str, Any]) -> None:
    try:
        with connect() as db:
            if job["job_type"] != "import_and_analyze" or not job.get("import_id"):
                raise ValueError(f"Unsupported job type: {job['job_type']}")
            summary = process_import(
                db,
                int(job["import_id"]),
                progress=lambda value: update_progress(int(job["job_id"]), value),
            )
        with connect() as db:
            db.execute(
                """UPDATE jobs SET status='complete', progress=1, metrics_json=?, updated_at=?
                   WHERE job_id=?""",
                (json_text(summary), utcnow(), job["job_id"]),
            )
    except Exception as exc:
        with connect() as db:
            db.execute(
                """UPDATE jobs SET status='failed', error=?, updated_at=? WHERE job_id=?""",
                (f"{type(exc).__name__}: {exc}", utcnow(), job["job_id"]),
            )
            if job.get("import_id"):
                db.execute(
                    "UPDATE imports SET import_status='failed' WHERE import_id=?",
                    (job["import_id"],),
                )


def main() -> None:
    initialize()
    signal.signal(signal.SIGTERM, _stop)
    signal.signal(signal.SIGINT, _stop)
    while running:
        job = claim_job()
        if job is None:
            time.sleep(0.75)
            continue
        run_job(job)


if __name__ == "__main__":
    main()

