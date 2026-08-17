from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


PALETTE = [
    "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
    "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
    "#393b79", "#637939", "#8c6d31", "#843c39", "#7b4173",
    "#3182bd", "#31a354", "#756bb1", "#e6550d", "#969696",
]


def _save(fig, path):
    fig.savefig(path, dpi=180, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def make_figures(result, outdir, dataset_kind):
    out = Path(outdir); figs = out / "figures"; figs.mkdir(parents=True, exist_ok=True)
    c = result.cells; q = c[c["pass_qc"] == True].copy()
    x = "x_centroid" if "x_centroid" in q else "x_centroid_um"; y = "y_centroid" if "y_centroid" in q else "y_centroid_um"

    fig, ax = plt.subplots(1, 3, figsize=(12, 3.4))
    ax[0].hist(c["computed_transcripts"], bins=35, color=PALETTE[0]); ax[0].set(title="Transcripts per cell", xlabel="Q20 gene transcripts")
    ax[1].hist(c["computed_features"], bins=35, color=PALETTE[2]); ax[1].set(title="Features per cell", xlabel="Detected genes")
    ax[2].hist(c["control_fraction"], bins=35, color=PALETTE[3]); ax[2].set(title="Control burden", xlabel="Control / all transcripts")
    fig.suptitle("Cell-level QC distributions"); _save(fig, figs / "01_cell_qc.png")

    point_size = max(0.25, min(10, 120000 / max(len(c), 1)))
    plot_n = min(250_000, len(c)); cp = c.sample(plot_n, random_state=17) if plot_n < len(c) else c
    fig, ax = plt.subplots(figsize=(7.2, 5.4)); sc=ax.scatter(cp[x], cp[y], c=np.log1p(cp["computed_transcripts"]), s=point_size, cmap="viridis", linewidths=0, rasterized=True); ax.invert_yaxis(); ax.set_aspect("equal"); ax.set(title="Spatial transcript yield", xlabel="x (µm)", ylabel="y (µm)"); fig.colorbar(sc, ax=ax, label="log1p transcripts/cell"); _save(fig, figs / "02_spatial_qc.png")

    labels = sorted(q["cell_type"].dropna().unique()); colors={v:PALETTE[i%len(PALETTE)] for i,v in enumerate(labels)}
    fig, ax = plt.subplots(figsize=(7.2, 5.4))
    for label in labels:
        z=q[q["cell_type"]==label]; z=z.sample(min(len(z), 40000), random_state=17); ax.scatter(z[x],z[y],s=point_size,label=label,color=colors[label],linewidths=0,rasterized=True)
    ax.invert_yaxis(); ax.set_aspect("equal"); ax.set(title="Marker-guided cell identities", xlabel="x (µm)", ylabel="y (µm)"); ax.legend(frameon=False,bbox_to_anchor=(1.02,1),loc="upper left",fontsize=8); _save(fig, figs / "03_cell_types.png")

    fig, ax=plt.subplots(figsize=(7,6)); im=ax.imshow(result.neighborhood, cmap="RdBu_r", vmin=-max(3,np.nanpercentile(abs(result.neighborhood),95)), vmax=max(3,np.nanpercentile(abs(result.neighborhood),95))); ax.set_xticks(range(len(result.neighborhood)),result.neighborhood.columns,rotation=45,ha="right",fontsize=8); ax.set_yticks(range(len(result.neighborhood)),result.neighborhood.index,fontsize=8); ax.set_title("Neighborhood enrichment (permutation z-score)"); fig.colorbar(im,ax=ax,label="z-score"); _save(fig, figs / "04_neighborhood_enrichment.png")

    top=result.gene_stats.head(15).sort_values("morans_i"); fig,ax=plt.subplots(figsize=(7.4,5.2)); ax.barh(top["gene"],top["morans_i"],color=PALETTE[0]); ax.axvline(0,color="#444",lw=.8); ax.set(title="Spatially coherent genes",xlabel="Moran's I"); _save(fig,figs/"05_morans_i.png")

    fig,ax=plt.subplots(figsize=(7.2,5.4)); ax.scatter(q["pc1"],q["pc2"],c=[colors[v] for v in q["cell_type"]],s=point_size,linewidths=0,rasterized=True); ax.set(title="Expression-state embedding",xlabel="PC/SVD 1",ylabel="PC/SVD 2"); _save(fig,figs/"06_embedding.png")

    proximity_types = [v for v in ["Tumor_epithelial", "T_cell", "B_cell", "Myeloid", "Fibroblast", "Endothelial"] if v in labels]
    if proximity_types and "distance_to_myoepithelial_um" in q:
        distance_groups = [q.loc[q["cell_type"].eq(v), "distance_to_myoepithelial_um"].dropna().clip(upper=1000) for v in proximity_types]
        fig, ax = plt.subplots(figsize=(8.2, 5.4))
        bp = ax.boxplot(distance_groups, tick_labels=proximity_types, showfliers=False, patch_artist=True)
        for patch, color in zip(bp["boxes"], [colors[v] for v in proximity_types]): patch.set_facecolor(color); patch.set_alpha(.75)
        ax.set(title="Proximity to myoepithelial cells", ylabel="Nearest-cell distance (µm; clipped at 1,000)")
        ax.tick_params(axis="x", rotation=35); ax.text(.01, .98, "Computational proxy—not a pathology-defined boundary", transform=ax.transAxes, va="top", color="#D65252")
        _save(fig, figs / "07_myoepithelial_proximity.png")
        q.groupby("cell_type")["distance_to_myoepithelial_um"].agg(["count", "median", "mean"]).sort_values("median").to_csv(out / "myoepithelial_proximity_summary.csv")

    q["cell_type"].value_counts().rename_axis("cell_type").rename("cells").to_csv(out / "cell_type_counts.csv")

    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    ax[0].hist(q["annotation_score"], bins=40, color=PALETTE[4]); ax[0].set(title="Annotation score", xlabel="Top marker-module score")
    ax[1].hist(q["annotation_margin"], bins=40, color=PALETTE[1]); ax[1].set(title="Annotation specificity", xlabel="Top minus second score")
    fig.suptitle("Marker annotation confidence"); _save(fig, figs / "08_annotation_confidence.png")

    counts = q["cell_type"].value_counts().sort_values()
    fig, ax = plt.subplots(figsize=(8, 5.5)); ax.barh(counts.index, counts.values, color=[colors[v] for v in counts.index]); ax.set(title="Cell-type composition", xlabel="QC-passing cells")
    _save(fig, figs / "09_cell_type_composition.png")

    if {"cell_area", "nucleus_cell_area_ratio", "transcripts_per_um2"}.issubset(q.columns):
        qp = q.sample(min(len(q), 100_000), random_state=17)
        fig, ax = plt.subplots(1, 2, figsize=(11, 4.5))
        ax[0].scatter(qp["cell_area"], qp["computed_transcripts"], s=1, alpha=.15, color=PALETTE[0], rasterized=True)
        ax[0].set(xscale="log", yscale="log", title="Yield versus segmented area", xlabel="Cell area (µm²)", ylabel="Gene transcripts")
        ax[1].hist(q["nucleus_cell_area_ratio"].dropna(), bins=50, color=PALETTE[2])
        ax[1].set(title="Nucleus-to-cell area ratio", xlabel="Nucleus area / cell area")
        fig.suptitle("Segmentation plausibility diagnostics"); _save(fig, figs / "10_segmentation_qc.png")

    pd.crosstab(q["cluster"], q["cell_type"], normalize="index").to_csv(out / "cluster_cell_type_fraction.csv")
    q.groupby("cell_type")[["annotation_score", "annotation_margin"]].agg(["count", "median", "mean"]).to_csv(out / "annotation_confidence_by_type.csv")

    result.cells.to_csv(out/"cells_with_qc_and_labels.csv.gz",compression="gzip")
    result.gene_stats.to_csv(out/"spatial_gene_statistics.csv",index=False)
    result.neighborhood.to_csv(out/"neighborhood_enrichment_z.csv")
    with (out/"qc_summary.json").open("w") as f: json.dump({**result.qc_summary,"dataset_kind":dataset_kind},f,indent=2)
