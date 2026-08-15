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

    fig, ax = plt.subplots(figsize=(7.2, 5.4)); sc=ax.scatter(c[x], c[y], c=np.log1p(c["computed_transcripts"]), s=10, cmap="viridis", linewidths=0); ax.invert_yaxis(); ax.set_aspect("equal"); ax.set(title="Spatial transcript yield", xlabel="x (µm)", ylabel="y (µm)"); fig.colorbar(sc, ax=ax, label="log1p transcripts/cell"); _save(fig, figs / "02_spatial_qc.png")

    labels = sorted(q["cell_type"].dropna().unique()); colors={v:PALETTE[i%len(PALETTE)] for i,v in enumerate(labels)}
    fig, ax = plt.subplots(figsize=(7.2, 5.4))
    for label in labels:
        z=q[q["cell_type"]==label]; ax.scatter(z[x],z[y],s=12,label=label,color=colors[label],linewidths=0)
    ax.invert_yaxis(); ax.set_aspect("equal"); ax.set(title="Marker-guided cell identities", xlabel="x (µm)", ylabel="y (µm)"); ax.legend(frameon=False,bbox_to_anchor=(1.02,1),loc="upper left",fontsize=8); _save(fig, figs / "03_cell_types.png")

    fig, ax=plt.subplots(figsize=(7,6)); im=ax.imshow(result.neighborhood, cmap="RdBu_r", vmin=-max(3,np.nanpercentile(abs(result.neighborhood),95)), vmax=max(3,np.nanpercentile(abs(result.neighborhood),95))); ax.set_xticks(range(len(result.neighborhood)),result.neighborhood.columns,rotation=45,ha="right",fontsize=8); ax.set_yticks(range(len(result.neighborhood)),result.neighborhood.index,fontsize=8); ax.set_title("Neighborhood enrichment (permutation z-score)"); fig.colorbar(im,ax=ax,label="z-score"); _save(fig, figs / "04_neighborhood_enrichment.png")

    top=result.gene_stats.head(15).sort_values("morans_i"); fig,ax=plt.subplots(figsize=(7.4,5.2)); ax.barh(top["gene"],top["morans_i"],color=PALETTE[0]); ax.axvline(0,color="#444",lw=.8); ax.set(title="Spatially coherent genes",xlabel="Moran's I"); _save(fig,figs/"05_morans_i.png")

    fig,ax=plt.subplots(figsize=(7.2,5.4)); ax.scatter(q["pc1"],q["pc2"],c=[colors[v] for v in q["cell_type"]],s=12,linewidths=0); ax.set(title="Expression-state embedding",xlabel="PC/SVD 1",ylabel="PC/SVD 2"); _save(fig,figs/"06_embedding.png")

    result.cells.to_csv(out/"cells_with_qc_and_labels.csv.gz",compression="gzip")
    result.gene_stats.to_csv(out/"spatial_gene_statistics.csv",index=False)
    result.neighborhood.to_csv(out/"neighborhood_enrichment_z.csv")
    with (out/"qc_summary.json").open("w") as f: json.dump({**result.qc_summary,"dataset_kind":dataset_kind},f,indent=2)
