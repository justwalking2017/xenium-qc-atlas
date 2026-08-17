from __future__ import annotations

import itertools
from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.cluster import MiniBatchKMeans
from sklearn.decomposition import TruncatedSVD
from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import normalize


CONTROL_PREFIXES = ("NegControl", "NegativeControl", "BLANK", "Unassigned")
NEGATIVE_FEATURE_TYPES = frozenset({"Negative Control Probe", "Negative Control Codeword"})


@dataclass
class AnalysisResult:
    cells: pd.DataFrame
    embedding: np.ndarray
    graph: sparse.csr_matrix
    gene_stats: pd.DataFrame
    neighborhood: pd.DataFrame
    qc_summary: dict
    marker_scores: pd.DataFrame


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
    # Prime 5K matrices also contain deprecated, unassigned and genomic-control
    # codewords. Only the two explicitly negative feature types estimate false
    # detection; folding all non-gene features together grossly over-filters.
    cmask = np.array([t in NEGATIVE_FEATURE_TYPES for t in feature_types])
    genomic_mask = np.array([t == "Genomic Control" for t in feature_types])
    X = matrix[:, gmask].astype(float)
    gene_names = np.array(genes)[gmask]
    totals = np.asarray(X.sum(axis=1)).ravel()
    detected = np.asarray((X > 0).sum(axis=1)).ravel()
    controls = np.asarray(matrix[:, cmask].sum(axis=1)).ravel() if cmask.any() else np.zeros(matrix.shape[0])
    genomic = np.asarray(matrix[:, genomic_mask].sum(axis=1)).ravel() if genomic_mask.any() else np.zeros(matrix.shape[0])
    cells["computed_transcripts"] = totals
    cells["computed_features"] = detected
    cells["control_fraction"] = controls / np.maximum(totals + controls, 1)
    cells["genomic_control_fraction"] = genomic / np.maximum(totals + genomic, 1)
    if "cell_area" in cells.columns:
        cells["transcripts_per_um2"] = totals / cells["cell_area"].to_numpy().clip(min=1e-6)
    pass_qc = (totals >= cfg["min_transcripts"]) & (detected >= cfg["min_features"]) & (cells["control_fraction"].to_numpy() <= cfg["max_control_fraction"])
    cells["pass_qc"] = pass_qc

    # Area-ratio and density flags expose segmentation failure modes rather than hiding them.
    if {"cell_area", "nucleus_area"}.issubset(cells.columns):
        cells["nucleus_cell_area_ratio"] = cells["nucleus_area"] / cells["cell_area"].clip(lower=1e-6)
        lo, hi = cells.loc[pass_qc, "nucleus_cell_area_ratio"].quantile([0.01, 0.99])
        cells["segmentation_outlier"] = ~cells["nucleus_cell_area_ratio"].between(lo, hi)
        density = cells.loc[pass_qc, "transcripts_per_um2"]
        dlo, dhi = density.quantile([0.01, 0.99])
        cells["density_outlier"] = ~cells["transcripts_per_um2"].between(dlo, dhi)
    else:
        cells["segmentation_outlier"] = False
        cells["density_outlier"] = False

    Xq = X[pass_qc]
    cq = cells.loc[pass_qc].copy()
    libnorm = normalize(Xq, norm="l1", axis=1) * 1e4
    logx = libnorm.copy(); logx.data = np.log1p(logx.data)
    n_comp = max(2, min(20, logx.shape[0] - 1, logx.shape[1] - 1))
    rng = np.random.default_rng(cfg["random_seed"])
    fit_n = min(int(cfg.get("max_fit_cells", 100_000)), logx.shape[0])
    fit_idx = rng.choice(logx.shape[0], fit_n, replace=False) if fit_n < logx.shape[0] else np.arange(logx.shape[0])
    svd = TruncatedSVD(n_components=n_comp, n_iter=5, random_state=cfg["random_seed"]).fit(logx[fit_idx])
    embedding = svd.transform(logx)
    n_clusters = max(2, min(cfg["n_clusters"], len(cq)))
    clusterer = MiniBatchKMeans(n_clusters=n_clusters, n_init=5, batch_size=4096, random_state=cfg["random_seed"])
    clusterer.fit(embedding[fit_idx])
    cq["cluster"] = clusterer.predict(embedding).astype(str)

    gene_to_idx = {g: i for i, g in enumerate(gene_names)}
    score_table = {}
    for label, marker_list in cfg["markers"].items():
        idx = [gene_to_idx[g] for g in marker_list if g in gene_to_idx]
        score_table[label] = np.asarray(logx[:, idx].mean(axis=1)).ravel() if idx else np.zeros(len(cq))
    scores = pd.DataFrame(score_table, index=cq.index)
    cq["cell_type"] = scores.idxmax(axis=1)
    sorted_scores = np.sort(scores.to_numpy(), axis=1)
    cq["annotation_score"] = sorted_scores[:, -1]
    cq["annotation_margin"] = sorted_scores[:, -1] - sorted_scores[:, -2]
    min_score = float(cfg.get("min_annotation_score", 0.05))
    min_margin = float(cfg.get("min_annotation_margin", 0.01))
    cq.loc[(cq["annotation_score"] < min_score) | (cq["annotation_margin"] < min_margin), "cell_type"] = "Unresolved"

    coords = cq[[x_col, y_col]].to_numpy(float)
    graph_n = min(int(cfg.get("max_graph_cells", 200_000)), len(cq))
    graph_idx = rng.choice(len(cq), graph_n, replace=False) if graph_n < len(cq) else np.arange(len(cq))
    graph_coords = coords[graph_idx]
    graph_cells = cq.iloc[graph_idx]
    k = max(1, min(cfg["k_neighbors"], len(cq) - 1))
    k = min(k, graph_n - 1)
    nbr = NearestNeighbors(n_neighbors=k + 1, n_jobs=-1).fit(graph_coords)
    indices = nbr.kneighbors(return_distance=False)[:, 1:]
    rows = np.repeat(np.arange(graph_n), k); cols = indices.ravel()
    graph = sparse.csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(graph_n, graph_n))
    graph = ((graph + graph.T) > 0).astype(float).tocsr()

    labels = sorted(cq["cell_type"].unique())
    label_to_code = {label: i for i, label in enumerate(labels)}
    codes = graph_cells["cell_type"].map(label_to_code).to_numpy(dtype=np.int32)
    n_labels = len(labels)
    source_codes = np.repeat(codes, k)
    target_index = indices.ravel()
    obs_array = np.bincount(
        source_codes * n_labels + codes[target_index], minlength=n_labels**2
    ).reshape(n_labels, n_labels)
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
    graph_logx = logx[graph_idx]
    means = np.asarray(graph_logx.mean(axis=0)).ravel(); sqmeans = np.asarray(graph_logx.power(2).mean(axis=0)).ravel(); variances = sqmeans - means**2
    for j in np.argsort(variances)[-min(200, len(gene_names)):]:
        v = graph_logx[:, j].toarray().ravel(); zc = v - v.mean(); denom = (zc @ zc)
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
        "density_outliers": int(cells["density_outlier"].sum()),
        "genes_in_panel": int(len(gene_names)),
        "cells_used_for_embedding_fit": int(fit_n),
        "cells_used_for_spatial_graph": int(graph_n),
        "unresolved_fraction": float(cq["cell_type"].eq("Unresolved").mean()),
    }
    # Preserve a focused set of vendor metrics as an orthogonal audit trail.
    if len(metrics):
        vendor_keys = (
            "fraction_transcripts_decoded_q20", "negative_control_probe_rate",
            "negative_control_codeword_rate", "adjusted_negative_control_probe_rate",
            "adjusted_negative_control_codeword_rate", "adjusted_genomic_control_probe_rate",
            "estimated_number_of_false_positive_transcripts_per_cell",
            "fraction_transcripts_assigned", "fraction_empty_cells",
            "segmented_cell_stain_frac", "segmented_cell_boundary_frac",
            "segmented_cell_interior_frac", "segmented_cell_nuc_expansion_frac",
        )
        qc_summary["xoa_metrics"] = {
            key: float(metrics.iloc[0][key]) for key in vendor_keys
            if key in metrics.columns and pd.notna(metrics.iloc[0][key])
        }
    return AnalysisResult(cells, embedding, graph, gene_stats, neighborhood, qc_summary, scores)
