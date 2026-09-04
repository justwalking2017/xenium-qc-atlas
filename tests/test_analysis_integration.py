import numpy as np
import pandas as pd
from scipy import sparse

from xenium_showcase.analysis import analyze


def test_analysis_emits_qc_anomaly_and_model_evaluation():
    rng = np.random.default_rng(3)
    n = 120
    barcodes = [f"cell-{i}" for i in range(n)]
    cells = pd.DataFrame({
        "cell_id": barcodes,
        "x_centroid": rng.uniform(0, 100, n),
        "y_centroid": rng.uniform(0, 100, n),
        "cell_area": rng.uniform(20, 80, n),
        "nucleus_area": rng.uniform(5, 18, n),
    })
    genes = ["EPCAM", "KRT8", "PTPRC", "CD3D", "NegControlProbe_1"]
    types = ["Gene Expression"] * 4 + ["Negative Control Probe"]
    values = rng.poisson(3, (n, len(genes)))
    matrix = sparse.csr_matrix(values)
    cfg = {
        "random_seed": 17, "min_transcripts": 2, "min_features": 1,
        "max_control_fraction": .5, "max_fit_cells": 100, "max_graph_cells": 100,
        "k_neighbors": 4, "n_clusters": 3, "n_permutations": 3,
        "min_annotation_score": .01, "min_annotation_margin": 0,
        "markers": {"Epithelial": ["EPCAM", "KRT8"], "T_cell": ["PTPRC", "CD3D"]},
        "qc": {"exclude_anomalies": False, "anomaly": {"enabled": True,
                "contamination": .05, "n_estimators": 20, "max_fit_cells": 100}},
    }
    result = analyze(cells, {"assignment_rate_all": .9, "assignment_rate_q20": .91,
                     "q20_fraction": .85, "transcripts_total": 2000, "median_qv": 40},
                     matrix, barcodes, genes, types, pd.DataFrame(), cfg)
    assert {"qc_rule_pass", "multivariate_anomaly", "anomaly_score", "qc_failure_reason"}.issubset(result.cells.columns)
    assert result.qc_summary["automatic_qc"]["method"] == "rules+isolation_forest"
    assert result.model_evaluation["scope"] == "internal_descriptive_evaluation"

