from __future__ import annotations

from copy import deepcopy


REQUIRED = {
    "project_name": str,
    "dataset_kind": str,
    "input_zip": str,
    "output_dir": str,
    "random_seed": int,
    "markers": dict,
}


def validate_config(cfg: dict) -> dict:
    """Validate and normalize a workflow configuration without hidden mutation."""
    if not isinstance(cfg, dict):
        raise ValueError("Configuration must be a YAML mapping")
    errors = []
    for key, expected in REQUIRED.items():
        if key not in cfg:
            errors.append(f"missing required key: {key}")
        elif not isinstance(cfg[key], expected):
            errors.append(f"{key} must be {expected.__name__}")
    qcfg = cfg.get("qc", {})
    if qcfg.get("mode", "hybrid") not in {"rules", "hybrid"}:
        errors.append("qc.mode must be one of: rules, hybrid")
    contamination = qcfg.get("anomaly", {}).get("contamination", 0.02)
    if not 0 < float(contamination) <= 0.5:
        errors.append("qc.anomaly.contamination must be in (0, 0.5]")
    if int(cfg.get("k_neighbors", 12)) < 1:
        errors.append("k_neighbors must be >= 1")
    if int(cfg.get("n_permutations", 50)) < 1:
        errors.append("n_permutations must be >= 1")
    if not cfg.get("markers"):
        errors.append("markers must contain at least one marker module")
    if errors:
        raise ValueError("Invalid workflow configuration:\n- " + "\n- ".join(errors))
    return deepcopy(cfg)
