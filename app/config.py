from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path


PARSER_VERSION = "backbone-v3"
STANDARDIZATION_VERSION = "lossless-v1"
MMP_ENGINE_VERSION = "strict-topology-v3"
RENDERER_VERSION = "rdkit-coordgen-svg-v5"


@dataclass(frozen=True)
class Settings:
    data_dir: Path
    database_path: Path
    asset_dir: Path

    @classmethod
    def from_env(cls) -> "Settings":
        data_dir = Path(os.getenv("PMMP_DATA_DIR", "data")).resolve()
        return cls(
            data_dir=data_dir,
            database_path=Path(os.getenv("PMMP_DB", data_dir / "pmmp.sqlite3")).resolve(),
            asset_dir=Path(os.getenv("PMMP_ASSET_DIR", data_dir / "assets")).resolve(),
        )


settings = Settings.from_env()

