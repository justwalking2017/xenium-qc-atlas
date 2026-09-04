# v0.3.0 — Automated QC and production workflow

- configurable hybrid rule-based QC and Isolation Forest anomaly screening;
- per-cell QC reason, anomaly flag and continuous anomaly score;
- versioned structured run metadata and JSON Schema;
- internal annotation/model diagnostics with explicit non-validation caveats;
- expanded unit and synthetic end-to-end tests plus Python-version CI matrix;
- Docker build and container smoke test in CI;
- complete assay, segmentation, model, spatial and inference failure-mode register;
- automatic anomaly map and report sections for QC, evaluation and limitations.

# v0.2.0 — Full breast-cancer biological run

This release upgrades the project from a technical smoke test to a complete biological Xenium analysis. It adds:

- CRC-validated analysis of the full official 10x FFPE human breast IDC output bundle;
- streaming QC over 86.8 million decoded transcripts;
- results for 574,852 segmented cells and 561,585 QC-passing cells;
- full-section cell identities, neighborhood enrichment, and Moran's I outputs;
- a myoepithelial-proximity application with an explicit pathology-proxy caveat;
- a new editable full-results PowerPoint deck;
- memory-efficient transcript parquet scanning and vectorized permutation counts.

The raw 22.7 GiB archive is not redistributed. Reproduction uses the official 10x URL and validates the expected byte length before analysis.

# v0.1.0 — Xenium QC-to-niche showcase

This release demonstrates practical Xenium analysis at three levels:

1. **Transcript layer:** Q-score-aware inputs, negative-control burden, transcript assignment, and spatial yield.
2. **Segmentation layer:** nucleus-to-cell area ratios, transcript density, segmentation outliers, and explicit retention of boundary-sensitive diagnostics.
3. **Biology layer:** marker-guided annotation, spatial kNN graphs, permutation-tested neighborhood enrichment, and Moran's I.

The bundled demo configuration targets the official 10x Xenium Prime mouse ileum tiny test dataset and is strictly a pipeline smoke test. The full breast-cancer configuration is intended for biological interpretation after downloading a complete public output bundle.

Reproducibility assets include pinned environment metadata, an MD5-verified download script, deterministic seeds, machine-readable QC outputs, figures, and an editable PowerPoint summary.
