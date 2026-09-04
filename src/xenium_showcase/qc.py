from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from sklearn.ensemble import IsolationForest
from sklearn.preprocessing import RobustScaler


@dataclass
class QCResult:
    pass_rules: np.ndarray
    anomaly: np.ndarray
    anomaly_score: np.ndarray
    reasons: pd.Series
    metadata: dict


def _reason_labels(failures: dict[str, np.ndarray], index) -> pd.Series:
    names = list(failures)
    values = np.full(len(index), "pass", dtype=object)
    for i in range(len(index)):
        failed = [name for name in names if failures[name][i]]
        if failed:
            values[i] = ";".join(failed)
    return pd.Series(values, index=index, name="qc_failure_reason")


def automatic_qc(cells: pd.DataFrame, cfg: dict, random_seed: int) -> QCResult:
    """Apply interpretable rules and optional multivariate anomaly detection.

    Rule failures determine ``pass_rules``. Anomaly calls are kept orthogonal by
    default, so a model cannot silently redefine the scientific analysis set.
    Set ``exclude_anomalies: true`` to include them in the final pass decision.
    """
    qcfg = cfg.get("qc", {})
    rules = qcfg.get("rules", {})
    totals = cells["computed_transcripts"].to_numpy(float)
    genes = cells["computed_features"].to_numpy(float)
    controls = cells["control_fraction"].to_numpy(float)
    failures = {
        "low_transcripts": totals < float(rules.get("min_transcripts", cfg.get("min_transcripts", 20))),
        "low_features": genes < float(rules.get("min_features", cfg.get("min_features", 10))),
        "high_control_fraction": controls > float(rules.get("max_control_fraction", cfg.get("max_control_fraction", 0.03))),
    }
    pass_rules = ~np.logical_or.reduce(list(failures.values()))

    acfg = qcfg.get("anomaly", {})
    anomaly = np.zeros(len(cells), dtype=bool)
    score = np.full(len(cells), np.nan)
    used_features: list[str] = []
    mode = str(qcfg.get("mode", "hybrid"))
    if mode == "hybrid" and acfg.get("enabled", True) and pass_rules.sum() >= 20:
        requested = acfg.get("features", [
            "computed_transcripts", "computed_features", "control_fraction",
            "cell_area", "nucleus_cell_area_ratio", "transcripts_per_um2",
        ])
        used_features = [name for name in requested if name in cells and cells[name].notna().any()]
        frame = cells.loc[pass_rules, used_features].replace([np.inf, -np.inf], np.nan)
        frame = frame.fillna(frame.median())
        positive = {"computed_transcripts", "computed_features", "cell_area", "transcripts_per_um2"}
        for name in positive.intersection(frame.columns):
            frame[name] = np.log1p(frame[name].clip(lower=0))
        values = RobustScaler().fit_transform(frame)
        rng = np.random.default_rng(random_seed)
        max_fit = min(int(acfg.get("max_fit_cells", 100_000)), len(values))
        fit_idx = rng.choice(len(values), max_fit, replace=False) if max_fit < len(values) else np.arange(len(values))
        model = IsolationForest(
            n_estimators=int(acfg.get("n_estimators", 200)),
            contamination=float(acfg.get("contamination", 0.02)),
            max_samples=min(int(acfg.get("max_samples", 10_000)), max_fit),
            random_state=random_seed,
            n_jobs=-1,
        ).fit(values[fit_idx])
        local_score = -model.score_samples(values)  # larger means more anomalous
        local_anomaly = model.predict(values) == -1
        positions = np.flatnonzero(pass_rules)
        score[positions] = local_score
        anomaly[positions] = local_anomaly

    reasons = _reason_labels(failures, cells.index)
    reasons.loc[anomaly & pass_rules] = "multivariate_anomaly"
    metadata = {
        "method": "rules+isolation_forest" if used_features else "rules_only",
        "mode": mode,
        "rule_failure_counts": {name: int(mask.sum()) for name, mask in failures.items()},
        "rule_pass_count": int(pass_rules.sum()),
        "anomaly_count_among_rule_pass": int((anomaly & pass_rules).sum()),
        "anomaly_features": used_features,
        "anomalies_excluded": bool(qcfg.get("exclude_anomalies", False)),
    }
    return QCResult(pass_rules, anomaly, score, reasons, metadata)
