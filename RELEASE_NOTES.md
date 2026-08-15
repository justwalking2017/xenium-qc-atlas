# v0.1.0 — Xenium QC-to-niche showcase

This release demonstrates practical Xenium analysis at three levels:

1. **Transcript layer:** Q-score-aware inputs, negative-control burden, transcript assignment, and spatial yield.
2. **Segmentation layer:** nucleus-to-cell area ratios, transcript density, segmentation outliers, and explicit retention of boundary-sensitive diagnostics.
3. **Biology layer:** marker-guided annotation, spatial kNN graphs, permutation-tested neighborhood enrichment, and Moran's I.

The bundled demo configuration targets the official 10x Xenium Prime mouse ileum tiny test dataset and is strictly a pipeline smoke test. The full breast-cancer configuration is intended for biological interpretation after downloading a complete public output bundle.

Reproducibility assets include pinned environment metadata, an MD5-verified download script, deterministic seeds, machine-readable QC outputs, figures, and an editable PowerPoint summary.

