# Workflow configuration

Each YAML profile is a complete, version-controlled analysis contract.

## Dataset and outputs

| Key | Meaning |
|---|---|
| `project_name` | Human-readable run name |
| `dataset_kind` | Intended role, such as biological showcase or format fixture |
| `dataset_url` | Public source/provenance URL |
| `input_zip` | Xenium output archive or extracted directory |
| `output_dir` | Generated artifacts; use a unique directory per run |
| `random_seed` | Controls all subsampling and fitted stochastic models |

## QC rules

`qc.rules` contains the interpretable inclusion thresholds. Legacy top-level
`min_transcripts`, `min_features`, and `max_control_fraction` remain supported.

```yaml
qc:
  mode: hybrid
  exclude_anomalies: false
  rules:
    min_transcripts: 20
    min_features: 10
    max_control_fraction: 0.03
```

Keep `exclude_anomalies: false` until anomaly maps have been reviewed against
histology and Xenium Explorer. If enabled, both rule failures and model calls
are excluded from embedding, annotation and spatial analysis.

Use `mode: rules` to disable model fitting while keeping the same transparent
rule engine. Use `mode: hybrid` to run both components.

## Anomaly model

```yaml
qc:
  anomaly:
    enabled: true
    method: isolation_forest
    contamination: 0.02
    n_estimators: 200
    max_samples: 10000
    max_fit_cells: 100000
    features: [computed_transcripts, computed_features, control_fraction,
               cell_area, nucleus_cell_area_ratio, transcripts_per_um2]
```

`contamination` sets the expected fraction of statistical outliers among cells
that pass hard rules; it is not an expected biological failure rate. Positive
skewed features are log-transformed, missing values are median-imputed, and all
features are robust-scaled before fitting.

## Expression, annotation and spatial analysis

| Key | Meaning |
|---|---|
| `max_fit_cells` | Seeded cells used to fit SVD and MiniBatchKMeans |
| `n_clusters` | Descriptive expression clusters; not cell types |
| `min_annotation_score` | Required top marker-module score |
| `min_annotation_margin` | Required top-minus-second module score |
| `markers` | Panel-aware marker modules |
| `max_graph_cells` | Seeded cells used in the physical spatial graph |
| `k_neighbors` | Nearest physical neighbors per graph node |
| `n_permutations` | Label permutations for neighborhood enrichment |

## Validation and reproducibility

```powershell
xenium-showcase --config configs/prime5k_human_breast.yml --validate-only
```

Every completed analysis writes `run_metadata.json`, `qc_summary.json`, and
`model_evaluation.json`. Archive hashes are recorded only for inputs at or below
the metadata hashing limit; very large public archives retain size, modification
time and source URL to avoid an expensive implicit read on every run.
