from __future__ import annotations

import numpy as np
import pandas as pd
from sklearn.metrics import normalized_mutual_info_score


def evaluate_models(cells: pd.DataFrame, marker_scores: pd.DataFrame, qc_metadata: dict) -> dict:
    """Return descriptive diagnostics; these are not external validation."""
    q = cells.loc[cells["pass_qc"] == True].copy()
    resolved = q["cell_type"].ne("Unresolved")
    cluster_nmi = (
        normalized_mutual_info_score(q.loc[resolved, "cluster"], q.loc[resolved, "cell_type"])
        if resolved.sum() > 1 else np.nan
    )
    per_type = {}
    for label, group in q.groupby("cell_type"):
        per_type[str(label)] = {
            "n": int(len(group)),
            "fraction": float(len(group) / max(len(q), 1)),
            "median_score": float(group["annotation_score"].median()),
            "median_margin": float(group["annotation_margin"].median()),
        }
    threshold_sensitivity = []
    for min_score in (0.0, 0.05, 0.10, 0.20):
        for min_margin in (0.0, 0.01, 0.05, 0.10):
            accepted = (q["annotation_score"] >= min_score) & (q["annotation_margin"] >= min_margin)
            threshold_sensitivity.append({
                "min_score": min_score, "min_margin": min_margin,
                "accepted_cells": int(accepted.sum()), "accepted_fraction": float(accepted.mean()),
            })
    return {
        "scope": "internal_descriptive_evaluation",
        "annotation": {
            "n_evaluated": int(len(q)),
            "resolved_fraction": float(resolved.mean()),
            "cluster_label_nmi_resolved": float(cluster_nmi),
            "score_quantiles": {str(k): float(v) for k, v in q["annotation_score"].quantile([.05, .25, .5, .75, .95]).items()},
            "margin_quantiles": {str(k): float(v) for k, v in q["annotation_margin"].quantile([.05, .25, .5, .75, .95]).items()},
            "per_type": per_type,
            "threshold_sensitivity": threshold_sensitivity,
        },
        "qc_anomaly": qc_metadata,
        "warnings": [
            "NMI measures cluster-label agreement, not biological correctness.",
            "No independent pathology labels or external reference atlas are available.",
            "Isolation Forest calls are unsupervised and sensitive to contamination and feature choice.",
        ],
    }
