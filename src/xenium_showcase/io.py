from __future__ import annotations

import gzip
import shutil
import zipfile
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from scipy import sparse


NEEDED = {"cells.parquet", "transcripts.parquet", "cell_feature_matrix.h5", "metrics_summary.csv"}


def materialize_input(source: str | Path, cache_dir: str | Path) -> Path:
    source, cache = Path(source), Path(cache_dir)
    if source.is_dir():
        return source
    if not source.exists():
        raise FileNotFoundError(f"Input not found: {source}")
    cache.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(source) as zf:
        names = {Path(n).name: n for n in zf.namelist()}
        missing = NEEDED - names.keys()
        if missing:
            raise ValueError(f"Xenium bundle lacks required files: {sorted(missing)}")
        for base in NEEDED:
            target = cache / base
            if not target.exists():
                with zf.open(names[base]) as src, target.open("wb") as dst:
                    shutil.copyfileobj(src, dst)
    return cache


def load_10x_h5(path: str | Path) -> tuple[sparse.csr_matrix, list[str], list[str], list[str]]:
    with h5py.File(path, "r") as h5:
        g = h5["matrix"]
        shape = tuple(g["shape"][:])
        matrix = sparse.csc_matrix((g["data"][:], g["indices"][:], g["indptr"][:]), shape=shape).T.tocsr()
        features = g["features"]
        names = [x.decode() for x in features["name"][:]]
        types = [x.decode() for x in features["feature_type"][:]] if "feature_type" in features else ["Gene Expression"] * len(names)
        barcodes = [x.decode() for x in g["barcodes"][:]]
    return matrix, barcodes, names, types


def summarize_transcripts(path: str | Path, batch_size: int = 1_000_000) -> dict:
    """Stream transcript QC without materializing tens of millions of string IDs."""
    pf = pq.ParquetFile(path)
    available = set(pf.schema_arrow.names)
    columns = [c for c in ("cell_id", "qv") if c in available]
    total = assigned = q20_total = q20_assigned = 0
    # QV is integer-like in Xenium output. A histogram avoids retaining tens of
    # millions of values merely to calculate a median.
    qv_hist = np.zeros(256, dtype=np.int64)
    for batch in pf.iter_batches(columns=columns, batch_size=batch_size):
        frame = batch.to_pandas()
        n = len(frame); total += n
        if "cell_id" in frame:
            is_assigned = frame["cell_id"].astype(str).ne("UNASSIGNED").to_numpy()
            assigned += int(is_assigned.sum())
        else:
            is_assigned = np.ones(n, dtype=bool)
        if "qv" in frame:
            qv = frame["qv"].to_numpy(dtype=np.float32, copy=False)
            high = qv >= 20
            q20_total += int(high.sum())
            q20_assigned += int((high & is_assigned).sum())
            bins = np.clip(np.rint(qv), 0, 255).astype(np.uint8)
            qv_hist += np.bincount(bins, minlength=256)
    if qv_hist.sum():
        median_qv = float(np.searchsorted(np.cumsum(qv_hist), (qv_hist.sum() + 1) // 2))
    else:
        median_qv = np.nan
    return {
        "transcripts_total": int(total),
        "assignment_rate_all": assigned / total if total else np.nan,
        "assignment_rate_q20": q20_assigned / q20_total if q20_total else np.nan,
        "q20_fraction": q20_total / total if total else np.nan,
        "median_qv": median_qv,
    }


def load_xenium(root: str | Path):
    root = Path(root)
    cells = pd.read_parquet(root / "cells.parquet")
    transcript_qc = summarize_transcripts(root / "transcripts.parquet")
    matrix, barcodes, genes, feature_types = load_10x_h5(root / "cell_feature_matrix.h5")
    metrics = pd.read_csv(root / "metrics_summary.csv")
    return cells, transcript_qc, matrix, barcodes, genes, feature_types, metrics
