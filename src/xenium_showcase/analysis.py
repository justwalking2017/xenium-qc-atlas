from __future__ import annotations

import itertools
from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.cluster import KMeans
from sklearn.decomposition import TruncatedSVD
from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import normalize


CONTROL_PREFIXES = ("NegControl", "NegativeControl", "BLANK", "Unassigned")


@dataclass
class AnalysisResult:
    cells: pd.DataFrame
    embedding: np.ndarray
    graph: sparse.csr_matrix
    gene_stats: pd.DataFrame
    neighborhood: pd.DataFrame
    qc_summary: dict


def _pick(df: pd.DataFrame, *candidates: str) -> str:
    for c in candidates:
        if c in df.columns:
            return c
    raise KeyError(f"None of {candidates} found in {list(df.columns)}")


def analyze(cells, transcript_qc, matrix, barcodes, genes, feature_types, metrics, cfg) -> AnalysisResult:
    cells = cells.copy()
    id_col = _pick(cells, "cell_id", "barcode")
    x_col = _pick(cells, "x_centroid", "x_centroid_um", "x")
    y_col = _pick(cells, "y_centroid", "y_centroid_um", "y")
    cells[id_col] = cells[id_col].astype(str)
    cells = cells.set_index(id_col).reindex(barcodes)

    gmask = np.array([t == "Gene Expression" and not n.startswith(CONTROL_PREFIXES) for n, t in zip(genes, feature_types)])
    cmask = ~gmask
    X = matrix[:, gmask].astype(float)
    gene_names = np.array(genes)[gmask]
    totals = np.asarray(X.sum(axis=1)).ravel()
    detected = np.asarray((X > 0).sum(axis=1)).ravel()
    controls = np.asarray(matrix[:, cmask].sum(axis=1)).ravel() if cmask.any() else np.zeros(matrix.shape[0])
    cells["computed_transcripts"] = totals
    cells["computed_features"] = detected
    cells["control_fraction"] = controls / np.maximum(totals + controls, 1)
    pass_qc = (totals >= cfg["min_transcripts"]) & (detected >= cfg["min_features"]) & (cells["control_fraction"].to_numpy() <= cfg["max_control_fraction"])
    cells["pass_qc"] = pass_qc

    # Area-ratio and density flags expose segmentation failure modes rather than hiding them.
    if {"cell_area", "nucleus_area"}.issubset(cells.columns):
        cells["nucleus_cell_area_ratio"] = cells["nucleus_area"] / cells["cell_area"].clip(lower=1e-6)
        lo, hi = cells.loc[pass_qc, "nucleus_cell_area_ratio"].quantile([0.01, 0.99])
        cells["segmentation_outlier"] = ~cells["nucleus_cell_area_ratio"].between(lo, hi)
    else:
        cells["segmentation_outlier"] = False

    Xq = X[pass_qc]
    cq = cells.loc[pass_qc].copy()
    libnorm = normalize(Xq, norm="l1", axis=1) * 1e4
    logx = libnorm.copy(); logx.data = np.log1p(logx.data)
    n_comp = max(2, min(20, logx.shape[0] - 1, logx.shape[1] - 1))
    embedding = TruncatedSVD(n_components=n_comp, random_state=cfg["random_seed"]).fit_transform(logx)
    n_clusters = max(2, min(cfg["n_clusters"], len(cq)))
    cq["cluster"] = KMeans(n_clusters=n_clusters, n_init=20, random_state=cfg["random_seed"]).fit_predict(embedding).astype(str)

    gene_to_idx = {g: i for i, g in enumerate(gene_names)}
    score_table = {}
    for label, marker_list in cfg["markers"].items():
        idx = [gene_to_idx[g] for g in marker_list if g in gene_to_idx]
        score_table[label] = np.asarray(logx[:, idx].mean(axis=1)).ravel() if idx else np.zeros(len(cq))
    scores = pd.DataFrame(score_table, index=cq.index)
    cq["cell_type"] = scores.idxmax(axis=1)
    cq["annotation_margin"] = np.sort(scores.to_numpy(), axis=1)[:, -1] - np.sort(scores.to_numpy(), axis=1)[:, -2]
    cq.loc[scores.max(axis=1) <= 0, "cell_type"] = "Unresolved"

    coords = cq[[x_col, y_col]].to_numpy(float)
    k = max(1, min(cfg["k_neighbors"], len(cq) - 1))
    nbr = NearestNeighbors(n_neighbors=k + 1).fit(coords)
    indices = nbr.kneighbors(return_distance=False)[:, 1:]
    rows = np.repeat(np.arange(len(cq)), k); cols = indices.ravel()
    graph = sparse.csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(len(cq), len(cq)))
    graph = ((graph + graph.T) > 0).astype(float).tocsr()

    labels = sorted(cq["cell_type"].unique())
    label_to_code = {label: i for i, label in enumerate(labels)}
    codes = cq["cell_type"].map(label_to_code).to_numpy(dtype=np.int32)
    n_labels = len(labels)
    source_codes = np.repeat(codes, k)
    target_index = indices.ravel()
    obs_array = np.bincount(
        source_codes * n_labels + codes[target_index], minlength=n_labels**2
    ).reshape(n_labels, n_labels)
    rng = np.random.default_rng(cfg["random_seed"])
    null = np.zeros((cfg["n_permutations"], len(labels), len(labels)))
    for p in range(cfg["n_permutations"]):
        shuffled = rng.permutation(codes)
        null[p] = np.bincount(
            np.repeat(shuffled, k) * n_labels + shuffled[target_index],
            minlength=n_labels**2,
        ).reshape(n_labels, n_labels)
    z = (obs_array - null.mean(0)) / np.maximum(null.std(0), 1)
    neighborhood = pd.DataFrame(z, index=labels, columns=labels)

    # Moran's I on expressed genes: an interpretable spatial-coherence check and application output.
    W = graph; wsum = W.sum(); gene_rows = []
    means = np.asarray(logx.mean(axis=0)).ravel(); sqmeans = np.asarray(logx.power(2).mean(axis=0)).ravel(); variances = sqmeans - means**2
    for j in np.argsort(variances)[-min(200, len(gene_names)):]:
        v = logx[:, j].toarray().ravel(); zc = v - v.mean(); denom = (zc @ zc)
        moran = (len(v) / wsum) * (zc @ (W @ zc)) / denom if denom > 0 and wsum else np.nan
        gene_rows.append((gene_names[j], means[j], variances[j], moran))
    gene_stats = pd.DataFrame(gene_rows, columns=["gene", "mean_log1p", "variance", "morans_i"]).sort_values("morans_i", ascending=False)

    cq["pc1"] = embedding[:, 0]; cq["pc2"] = embedding[:, 1]
    myo_mask = cq["cell_type"].eq("Myoepithelial").to_numpy()
    if myo_mask.any():
        myo_nn = NearestNeighbors(n_neighbors=1).fit(coords[myo_mask])
        cq["distance_to_myoepithelial_um"] = myo_nn.kneighbors(coords, return_distance=True)[0].ravel()
    else:
        cq["distance_to_myoepithelial_um"] = np.nan
    cells.loc[cq.index, cq.columns] = cq
    qc_summary = {
        "cells_total": int(len(cells)), "cells_pass_qc": int(pass_qc.sum()), "pass_rate": float(pass_qc.mean()),
        "median_transcripts": float(np.median(totals)), "median_features": float(np.median(detected)),
        "median_control_fraction": float(np.median(cells["control_fraction"])),
        "transcript_assignment_rate": transcript_qc["assignment_rate_all"],
        "q20_assignment_rate": transcript_qc["assignment_rate_q20"],
        "q20_transcript_fraction": transcript_qc["q20_fraction"],
        "transcripts_total": transcript_qc["transcripts_total"],
        "median_qv": transcript_qc["median_qv"],
        "segmentation_outliers": int(cells["segmentation_outlier"].sum()),
    }
    return AnalysisResult(cells, embedding, graph, gene_stats, neighborhood, qc_summary)
