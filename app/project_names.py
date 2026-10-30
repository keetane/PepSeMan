"""Consistent, collision-free names based on creation date and source filename."""

from __future__ import annotations

import re
import sqlite3
from datetime import datetime
from pathlib import PurePosixPath
from zoneinfo import ZoneInfo


PROJECT_TIMEZONE = ZoneInfo("Asia/Tokyo")


def project_name(filename: str, created_at: str) -> str:
    basename = PurePosixPath(filename.replace("\\", "/")).name
    stem = PurePosixPath(basename).stem or basename
    if not PurePosixPath(basename).suffix:
        stem = re.sub(r"^\d{8}_", "", stem)
    stem = re.sub(r"[^\w-]+", "_", stem, flags=re.UNICODE).strip("_-") or "project"
    date = datetime.fromisoformat(created_at).astimezone(PROJECT_TIMEZONE).strftime("%Y%m%d")
    return f"{date}_{stem[:111]}"


def unique_project_name(db: sqlite3.Connection, filename: str, created_at: str, *, exclude_id: int | None = None) -> str:
    base = project_name(filename, created_at)
    existing = {
        row["name"] for row in db.execute(
            "SELECT name FROM projects WHERE project_id != ?", (exclude_id or -1,)
        )
    }
    candidate = base
    suffix = 2
    while candidate in existing:
        tail = f"_{suffix}"
        candidate = f"{base[:120 - len(tail)]}{tail}"
        suffix += 1
    return candidate
