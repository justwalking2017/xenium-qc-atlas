from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


PALETTE = ["#3D8DFF", "#EF5DA8", "#15A66D", "#FF9D2E", "#7657D5", "#00A6A6", "#D65252", "#6B7280"]


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
    fig, ax = plt.subplots(figsize=(7.2, 5.4)); sc=ax.scatter(c[x], c[y], c=np.log1p(c["computed_transcripts"]), s=point_size, cmap="viridis", linewidths=0, rasterized=True); ax.invert_yaxis(); ax.set_aspect("equal"); ax.set(title="Spatial transcript yield", xlabel="x (µm)", ylabel="y (µm)"); fig.colorbar(sc, ax=ax, label="log1p transcripts/cell"); _save(fig, figs / "02_spatial_qc.png")

    labels = sorted(q["cell_type"].dropna().unique()); colors={v:PALETTE[i%len(PALETTE)] for i,v in enumerate(labels)}
    fig, ax = plt.subplots(figsize=(7.2, 5.4))
    for label in labels:
        z=q[q["cell_type"]==label]; ax.scatter(z[x],z[y],s=point_size,label=label,color=colors[label],linewidths=0,rasterized=True)
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

    result.cells.to_csv(out/"cells_with_qc_and_labels.csv.gz",compression="gzip")
    result.gene_stats.to_csv(out/"spatial_gene_statistics.csv",index=False)
    result.neighborhood.to_csv(out/"neighborhood_enrichment_z.csv")
    with (out/"qc_summary.json").open("w") as f: json.dump({**result.qc_summary,"dataset_kind":dataset_kind},f,indent=2)
