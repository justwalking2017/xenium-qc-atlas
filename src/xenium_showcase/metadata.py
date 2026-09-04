from __future__ import annotations

import hashlib
import json
import platform
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy
import pandas
import scipy
import sklearn


def _git_value(*args: str) -> str | None:
    try:
        return subprocess.check_output(["git", *args], text=True, stderr=subprocess.DEVNULL).strip()
    except (OSError, subprocess.CalledProcessError):
        return None


def file_identity(path: str | Path) -> dict:
    path = Path(path)
    info = {"path": str(path), "exists": path.exists()}
    if path.exists():
        info.update({"size_bytes": path.stat().st_size, "modified_utc": datetime.fromtimestamp(path.stat().st_mtime, timezone.utc).isoformat()})
        if path.is_file() and path.stat().st_size <= 100_000_000:
            digest = hashlib.sha256()
            with path.open("rb") as handle:
                for block in iter(lambda: handle.read(1024 * 1024), b""):
                    digest.update(block)
            info["sha256"] = digest.hexdigest()
    return info


def build_run_metadata(cfg: dict, config_path: str | Path, output_dir: str | Path) -> dict:
    return {
        "schema_version": "1.0.0",
        "run": {
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "project_name": cfg.get("project_name"),
            "dataset_kind": cfg.get("dataset_kind"),
            "random_seed": cfg.get("random_seed"),
            "output_dir": str(output_dir),
        },
        "dataset": {"url": cfg.get("dataset_url"), "input": file_identity(cfg.get("input_zip", ""))},
        "workflow": {"config": file_identity(config_path), "configuration": cfg},
        "software": {
            "python": sys.version.split()[0], "platform": platform.platform(),
            "numpy": numpy.__version__, "pandas": pandas.__version__, "scipy": scipy.__version__,
            "scikit_learn": sklearn.__version__, "git_commit": _git_value("rev-parse", "HEAD"),
            "git_dirty": bool(_git_value("status", "--porcelain")),
        },
    }


def write_json(path: str | Path, value: dict) -> None:
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False), encoding="utf-8")

