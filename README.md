# Xenium QC Atlas: from decoded transcripts to tissue niches

[![Python](https://img.shields.io/badge/python-3.10%2B-3776AB)](https://www.python.org/) [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE) [![Data: CC BY 4.0](https://img.shields.io/badge/demo%20data-CC%20BY%204.0-green.svg)](https://www.10xgenomics.com/support/software/xenium-onboard-analysis/latest/resources/xenium-example-data)

**Reviewer-ready deliverable:** [Open or download the self-contained Prime 5K QC + analysis report](Xenium_Prime5K_QC_Analysis_Report.html).

A release-ready, hands-on Xenium project that treats quality control as a chain of evidence: **decoding → transcript assignment → segmentation → cell profiles → spatial biology**. It produces auditable tables, publication-ready figures, and a self-contained HTML report.

## Production engineering features

- **Configurable workflow:** YAML profiles control QC rules, anomaly features, sampling, clustering, annotation and spatial permutations; `--validate-only` checks the profile and input before a long run.
- **Automatic QC:** transparent transcript, feature and control rules emit per-cell failure reasons.
- **Anomaly detection:** a seeded Isolation Forest screens joint yield, control, geometry and density abnormalities. Calls are advisory unless `qc.exclude_anomalies: true` is explicit.
- **Tests and CI:** unit and synthetic end-to-end tests run on Python 3.10–3.12; CI also builds and smoke-tests the container.
- **Containerization:** code and dependencies are packaged while public/raw data are mounted at runtime.
- **Structured metadata:** each run records dataset identity, resolved configuration, software versions, git commit and dirty-tree state under a versioned JSON contract.
- **Model evaluation:** annotation score/margin distributions, resolved fraction and cluster-label NMI are separated from biological results.
- **Failure-mode register:** assay, segmentation, annotation, anomaly-model, spatial-statistics and study-design risks include symptoms and mitigations.

## Why this is more than a standard single-cell workflow

Xenium is targeted and image-based. A credible analysis therefore cannot rely only on total counts and genes per cell. This project checks:

- assay specificity using negative-control features and transcript QV;
- transcript-to-cell assignment and local transcript yield;
- segmentation plausibility using cell/nucleus geometry and boundary-sensitive outliers;
- panel-aware annotation confidence instead of pretending a targeted panel is whole transcriptome;
- spatial coherence and cell–cell adjacency against permutation nulls;
- failure modes by retaining unresolved cells and reporting sensitivity to thresholds.

The paper suggested in the project brief, **Moses et al., Nature Methods (2024)**, benchmarks sequencing-based spatial methods. Its principles—orthogonal QC axes, sensitivity/specificity, spatial concordance, and reproducibility—are useful, but the implementation here adapts them to Xenium-specific decoded transcripts and cell segmentation. Platform-specific choices are further grounded in Janesick et al. (2023) and Goods et al. (2025).

## Analysis profiles

| Profile | Role | Dataset / panel | Biological interpretation |
|---|---|---|---|
| `prime5k_human_breast.yml` | **Primary reviewer-ready analysis** | Independent FFPE breast cancer specimen; Prime 5K Human Pan Tissue & Pathways + 100 custom genes | Main QC, annotation, neighborhood, Moran's I, and segmentation showcase |
| `full_human_breast.yml` | Legacy Xenium v1 comparison | Independent FFPE breast IDC specimen; smaller Xenium v1 breast panel + add-on genes | Historical workflow and panel-depth comparison; not a replicate of the Prime 5K specimen |
| `demo_mouse_ileum.yml` | Format-test fixture only | Artificially trimmed Prime 5K mouse ileum patches (~12 MB) | **No biological claims**; validates XOA file compatibility and the code path |

> **Important:** the two breast profiles are independent public 10x datasets generated with different panels and XOA versions. Their cell counts and QC metrics are not repeated measurements of the same specimen and should not be compared as replicates.

## Quick start

### Primary workflow: Prime 5K breast cancer

Create the environment, download the complete public dataset, and run:

```powershell
conda env create -f environment.yml
conda activate xenium-showcase
.\scripts\download_prime5k_breast.ps1
xenium-showcase --config configs\prime5k_human_breast.yml
```

Validate the workflow and archive without running the expensive analysis:

```powershell
xenium-showcase --config configs\prime5k_human_breast.yml --validate-only
```

This profile downloads the official 38.17 GiB output bundle, selectively extracts
the files required for computation, and writes a self-contained QC + analysis
report to `results/human_breast_prime5k/xenium_prime5k_qc_analysis_report.html`.

Expected outputs:

```text
results/human_breast_prime5k/
├── figures/01_cell_qc.png ... 10_segmentation_qc.png
├── qc_summary.json
├── model_evaluation.json
├── run_metadata.json
├── cells_with_qc_and_labels.csv.gz
├── spatial_gene_statistics.csv
├── neighborhood_enrichment_z.csv
└── xenium_prime5k_qc_analysis_report.html
```

### Container

```powershell
docker build -t xenium-qc-atlas:0.3.0 .
```

Mount data and result directories at runtime and provide a config whose input/output paths match those mount points. Raw data, generated results, PDFs and slide decks are deliberately excluded from the image.

Read the [complete limitations and failure-mode register](docs/limitations_and_failure_modes.md) before interpreting anomaly, annotation or spatial-statistics outputs.

### Optional format test

Use the tiny mouse ileum fixture only to verify installation and file compatibility:

```powershell
powershell -ExecutionPolicy Bypass -File scripts/download_demo.ps1
xenium-showcase --config configs/demo_mouse_ileum.yml
```

## Primary executed result: Prime 5K breast cancer

The complete 38.17 GiB 10x FFPE human breast cancer output bundle was downloaded, byte-verified, selectively extracted, and analyzed end to end:

| Metric | Prime 5K result |
|---|---:|
| XOA cell objects | 699,110 |
| Cells passing configured QC | 504,483 (72.2%) |
| Gene-expression features | 5,101 |
| Transcripts in transcript parquet | 109,411,890 |
| Q20 transcript fraction | 85.2% |
| Q20 transcript assignment | 89.0% |
| Median gene transcripts / cell | 49 |
| Median detected genes / cell | 46 |
| Conservatively unresolved passing cells | 179,757 (35.6%) |

The 72.2% pass rate is driven almost entirely by cell objects with fewer than 20 gene transcripts, not by elevated negative-control signal. `Unresolved` is an annotation-confidence outcome rather than a QC failure: most unresolved cells lack signal from the deliberately compact marker dictionary, while a smaller group has ambiguous top-two module scores.

Spatially coherent candidate genes included `CA12`, `XBP1`, `RAB11FIP1`, `TSPAN13`, `ANKRD30A`, `CCND1`, `GATA3`, and `ESR1`. Neighborhood enrichment and Moran's I are treated as hypothesis-generating outputs that require pathology regions and patient-level validation.

## Legacy comparison: Xenium v1 breast panel

Before the Prime 5K upgrade, the workflow was executed on a separate 10x FFPE human breast IDC dataset using the smaller Xenium v1 breast panel. These metrics are retained to document project evolution and provide a panel-depth comparison:

| Metric | Result |
|---|---:|
| Segmented cells | 574,852 |
| Cells passing configured QC | 561,585 (97.7%) |
| Decoded transcripts in transcript parquet | 86,844,040 |
| Q20 transcript fraction | 78.5% |
| Median transcripts / cell | 105 |
| Median detected genes / cell | 46 |
| Segmentation geometry outliers | 11,562 |

Marker-guided annotation recovered 316,693 tumor epithelial, 133,590 fibroblast, 33,592 myoepithelial, 33,013 T, 24,265 endothelial, 9,306 myeloid, and 7,793 B cells. Spatially coherent genes included `CLIC6`, `SERPINA3`, `GATA3`, `PGR`, `KRT14`, `LUM`, and `POSTN`.

The myoepithelial-proximity analysis is explicitly treated as a computational nearest-cell proxy. It is not a substitute for pathology-defined DCIS or invasive boundaries.

Run the legacy v1 profile after downloading its independent official bundle:

```powershell
powershell -ExecutionPolicy Bypass -File scripts/download_full_breast.ps1
xenium-showcase --config configs/full_human_breast.yml
```

## QC decision framework

| Layer | Metric | Interpretation | Action |
|---|---|---|---|
| Decoding | transcript QV, negative controls | specificity / false calls | inspect tails and spatial hotspots before filtering |
| Assignment | assigned fraction, transcripts/cell | boundary performance | compare by tissue region and cell morphology |
| Segmentation | nucleus:cell area, transcript density | over/under-segmentation | flag extremes; review in Xenium Explorer |
| Cell profile | counts, features, annotation margin | panel coverage and identity | use cell-type-aware thresholds; retain unresolved |
| Spatial | local yield, Moran's I, adjacency z | artifacts vs coherent biology | require morphology or replicate support |

Thresholds are configuration values, not universal truths. For real cohorts, derive them per sample and cell class, then perform sensitivity analysis. Avoid deleting low-RNA immune cells or large adipocytes with a single global cutoff.

## Application story: breast cancer boundary ecology

The primary Prime 5K profile is designed around a question Xenium is unusually good at answering: **how myoepithelial integrity and immune/stromal neighborhoods change across DCIS-to-invasive boundaries**. The analysis can be extended with:

1. pathology-defined DCIS, invasive, and normal-adjacent regions;
2. distance-to-boundary gradients for ACTA2/KRT15 myoepithelial cells;
3. permutation-tested enrichment of tumor–myeloid, tumor–fibroblast, and B–T neighborhoods;
4. spatially coherent marker programs and region-aware pseudobulk differential testing;
5. sensitivity analyses across segmentation versions or nucleus expansion radii.

This direction follows recent Xenium applications emphasizing treatment-responsive immune niches, residual immune-cold tumor pockets, and spatially co-infiltrating fibroblast–macrophage programs—not just cluster maps.

## What I would add for a production cohort

- replicate-aware pseudobulk models (patient is the unit of inference);
- pathology registration and blinded ROI definitions;
- segmentation benchmarking against membrane staining or alternative algorithms;
- ambient/spillover diagnostics and mixed-lineage scores near dense boundaries;
- versioned panel, XOA, segmentation, and coordinate-system provenance;
- Nextflow/Snakemake orchestration plus container digests for multi-sample scale.

## References

- Moses L. et al. Systematic comparison of sequencing-based spatial transcriptomic methods. *Nature Methods* (2024). https://doi.org/10.1038/s41592-024-02325-3
- Janesick A. et al. High resolution mapping of the tumor microenvironment using integrated single-cell, spatial and in situ analysis. *Nature Communications* (2023). https://doi.org/10.1038/s41467-023-43458-x
- Goods B.A. et al. Optimizing Xenium In Situ data utility by quality assessment and best-practice analysis workflows. *Nature Methods* (2025). https://doi.org/10.1038/s41592-025-02617-2
- Sakai S.A. et al. Single-cell spatial analysis with Xenium reveals anti-tumour responses after radiotherapy plus anti-PD-L1. *British Journal of Cancer* (2025). https://doi.org/10.1038/s41416-025-03088-0

## License and data provenance

Code is MIT licensed. 10x public datasets remain under their stated CC BY 4.0 terms and are not redistributed in this repository. Cite 10x Genomics and the source dataset when publishing figures.
