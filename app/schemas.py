from __future__ import annotations

from typing import Any, Literal

from pydantic import BaseModel, Field, field_validator


class ProjectCreate(BaseModel):
    name: str = Field(min_length=1, max_length=120)

    @field_validator("name")
    @classmethod
    def strip_name(cls, value: str) -> str:
        return value.strip()


class ProjectUpdate(ProjectCreate):
    pass


class Project(BaseModel):
    project_id: int
    name: str
    settings_json: str
    created_at: str
    updated_at: str


class ImportResponse(BaseModel):
    import_id: int
    job_id: int | None = None
    duplicate: bool = False
    dataset_revision: int


class SortSpec(BaseModel):
    column: str
    direction: Literal["asc", "desc"] = "asc"


class Page(BaseModel):
    rows: list[dict[str, Any]]
    total: int
    page: int
    page_size: int

