# Limitations and failure modes

This register is part of the analysis contract. A successful run is not evidence that every biological conclusion is valid.

## Assay and panel

| Failure mode | Observable symptom | Consequence | Detection | Mitigation |
|---|---|---|---|---|
| Targeted-panel blind spot | Expected identity has few/no panel markers | Cells become unresolved or mislabelled | Marker coverage table; unresolved rate | Add probes; validate with an external assay/reference |
| Probe-specific background | Signal correlates with negative/genomic controls | False expression and spatial patterns | Control maps and feature-type audit | Exclude affected probes; review XOA metrics |
| Low RNA preservation | Broad low counts/features | Loss of sensitivity and selective cell dropout | Spatial yield maps; histology concordance | Optimize fixation/preanalytics; stratify interpretation |
| Optical crowding/decode saturation | High-density regions plateau or lose QV | Biased abundance estimates | QV and density maps; XOA saturation metrics | Treat counts as semi-quantitative; orthogonal validation |

## Segmentation and assignment

| Failure mode | Observable symptom | Consequence | Detection | Mitigation |
|---|---|---|---|---|
| Over-segmentation | Tiny cells, low counts, extreme nucleus/cell ratio | One cell split into fragments | Area, density and geometry outliers | Re-segment; merge only with validated rules |
| Under-segmentation | Very large/high-count objects | Multiple cells merged | Area–yield plot; morphology | Re-segment; do not interpret as a single cell |
| Boundary expansion contamination | Low density or mixed marker modules | Neighbor transcripts assigned to a cell | Density, margin, morphology | Adjust expansion; transcript-level analysis |
| Unassigned transcripts | Spatially concentrated extracellular signal | Missing cell expression | Q20 assignment map/rate | Review segmentation and tissue morphology |

## QC anomaly model

Isolation Forest is unsupervised. It learns rarity, not biological invalidity. Rare, valid cell types may be called anomalous; common systematic artifacts may appear normal. Its output is therefore advisory by default. Results depend on contamination, feature transformations, sampling and dataset composition. Refit per dataset and inspect anomaly maps before exclusion.

## Annotation

Marker-module annotation is not ground truth. Scores are panel-dependent and do not model ambient/neighbor contamination. The highest score can still be wrong. `Unresolved` combines low information, mixed states, missing marker coverage and true intermediate biology. Cluster-label NMI measures internal agreement only; high NMI does not prove biological accuracy.

## Spatial statistics

KNN results depend on `k`, sampling density and tissue geometry. Label permutation tests spatial randomness conditional on observed labels and fixed coordinates; it does not test causal interaction. Moran's I ranks autocorrelation but the current gene table has no gene-wise permutation p-values or FDR. Spatially adjacent cells are not independent replicates.

## Study design and inference

One section from one patient supports within-section exploration only. Cells cannot substitute for patients. Multi-section studies require section-aware processing and patient-level or hierarchical inference. Batch, ROI selection, sectioning depth and tissue composition can confound comparisons. Clinical claims require independent cohorts and pathology/orthogonal validation.

## Operational failures

- Interrupted or partial archive extraction: run `--validate-only`; required members are checked.
- Barcode mismatch: cells are reindexed to matrix barcodes; missing IDs must be treated as an input integrity failure.
- Memory exhaustion: lower `max_fit_cells` and `max_graph_cells`; this changes computational precision and must be recorded.
- Non-reproducible sampling: preserve `random_seed`, config, software versions and git commit in `run_metadata.json`.
- HTML opened from an incomplete result directory: use the self-contained report generated at the end of a successful run.

