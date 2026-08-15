from __future__ import annotations

import gzip
import shutil
import zipfile
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy import sparse


NEEDED = {"cells.parquet", "transcripts.parquet", "cell_feature_matrix.h5", "metrics_summary.csv", "gene_panel.json"}


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


def load_xenium(root: str | Path):
    root = Path(root)
    cells = pd.read_parquet(root / "cells.parquet")
    transcripts = pd.read_parquet(root / "transcripts.parquet")
    matrix, barcodes, genes, feature_types = load_10x_h5(root / "cell_feature_matrix.h5")
    metrics = pd.read_csv(root / "metrics_summary.csv")
    return cells, transcripts, matrix, barcodes, genes, feature_types, metrics

