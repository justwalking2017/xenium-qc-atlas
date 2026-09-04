import numpy as np
import pandas as pd

from xenium_showcase.qc import automatic_qc


def _cells(n=100):
    rng = np.random.default_rng(7)
    frame = pd.DataFrame({
        "computed_transcripts": rng.poisson(100, n).astype(float),
        "computed_features": rng.poisson(50, n).astype(float),
        "control_fraction": rng.uniform(0, .005, n),
        "cell_area": rng.lognormal(3.5, .2, n),
        "nucleus_cell_area_ratio": rng.normal(.45, .05, n),
    })
    frame["transcripts_per_um2"] = frame["computed_transcripts"] / frame["cell_area"]
    return frame


def test_rules_are_auditable_and_non_overlapping_with_model_decision():
    cells = _cells()
    cells.loc[0, "computed_transcripts"] = 1
    cells.loc[1, "computed_features"] = 1
    cells.loc[2, "control_fraction"] = .5
    cfg = {"min_transcripts": 20, "min_features": 10, "max_control_fraction": .03,
           "qc": {"exclude_anomalies": False, "anomaly": {"enabled": False}}}
    result = automatic_qc(cells, cfg, 17)
    assert result.pass_rules.sum() == 97
    assert result.reasons.iloc[0] == "low_transcripts"
    assert result.reasons.iloc[1] == "low_features"
    assert result.reasons.iloc[2] == "high_control_fraction"
    assert not result.anomaly.any()


def test_isolation_forest_is_reproducible():
    cells = _cells(250)
    cfg = {"qc": {"anomaly": {"enabled": True, "contamination": .04,
                                "n_estimators": 25, "max_fit_cells": 200}}}
    first = automatic_qc(cells, cfg, 17)
    second = automatic_qc(cells, cfg, 17)
    np.testing.assert_array_equal(first.anomaly, second.anomaly)
    np.testing.assert_allclose(first.anomaly_score, second.anomaly_score, equal_nan=True)

