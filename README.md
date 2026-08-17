# Xenium QC Atlas: from decoded transcripts to tissue niches

[![Python](https://img.shields.io/badge/python-3.10%2B-3776AB)](https://www.python.org/) [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE) [![Data: CC BY 4.0](https://img.shields.io/badge/demo%20data-CC%20BY%204.0-green.svg)](https://www.10xgenomics.com/support/software/xenium-onboard-analysis/latest/resources/xenium-example-data)

**Reviewer-ready deliverable:** [Open or download the self-contained Prime 5K QC + analysis report](Xenium_Prime5K_QC_Analysis_Report.html).

A release-ready, hands-on Xenium project that treats quality control as a chain of evidence: **decoding → transcript assignment → segmentation → cell profiles → spatial biology**. It produces auditable tables, publication-ready figures, and an editable presentation.

## Why this is more than a standard single-cell workflow

Xenium is targeted and image-based. A credible analysis therefore cannot rely only on total counts and genes per cell. This project checks:

- assay specificity using negative-control features and transcript QV;
- transcript-to-cell assignment and local transcript yield;
- segmentation plausibility using cell/nucleus geometry and boundary-sensitive outliers;
- panel-aware annotation confidence instead of pretending a targeted panel is whole transcriptome;
- spatial coherence and cell–cell adjacency against permutation nulls;
- failure modes by retaining unresolved cells and reporting sensitivity to thresholds.

The paper suggested in the project brief, **Moses et al., Nature Methods (2024)**, benchmarks sequencing-based spatial methods. Its principles—orthogonal QC axes, sensitivity/specificity, spatial concordance, and reproducibility—are useful, but the implementation here adapts them to Xenium-specific decoded transcripts and cell segmentation. Platform-specific choices are further grounded in Janesick et al. (2023) and Goods et al. (2025).

## Two execution profiles

| Profile | Purpose | Size | Biological claims |
|---|---|---:|---|
| `demo_mouse_ileum.yml` | Official 10x trimmed output; tests XOA formats and the complete code path | ~12 MB | **No** — only a few artificial patches |
| `full_human_breast.yml` | Complete FFPE breast section; tumor–myoepithelial boundary and immune niches | tens of GB | Yes, with pathology/replicate validation |

## Quick start

For the complete Prime 5K breast-cancer workflow, download and run:

```powershell
.\scripts\download_prime5k_breast.ps1
conda run -n xenium-qc-atlas xenium-showcase --config configs\prime5k_human_breast.yml
```

This profile downloads the official 38.17 GiB output bundle, selectively extracts
the files required for computation, and writes a self-contained QC + analysis
report to `results/human_breast_prime5k/xenium_prime5k_qc_analysis_report.html`.

```powershell
conda env create -f environment.yml
conda activate xenium-showcase
powershell -ExecutionPolicy Bypass -File scripts/download_demo.ps1
xenium-showcase --config configs/demo_mouse_ileum.yml
```

Expected outputs:

```text
results/<profile>/
├── figures/01_cell_qc.png ... 06_embedding.png
├── qc_summary.json
├── cells_with_qc_and_labels.csv.gz
├── spatial_gene_statistics.csv
└── neighborhood_enrichment_z.csv
```

## Executed full breast-cancer run

The complete 10x FFPE human breast IDC output bundle was downloaded, CRC-validated, and analyzed end to end. The checked-in results were generated from the full biological sample, not the trimmed format fixture:

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

Run the same full profile after downloading the official bundle:

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

The full profile is designed around a question Xenium is unusually good at answering: **how myoepithelial integrity and immune/stromal neighborhoods change across DCIS-to-invasive boundaries**. The analysis can be extended with:

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
