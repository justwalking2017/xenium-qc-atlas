from __future__ import annotations

import argparse
from pathlib import Path

import yaml

from .analysis import analyze
from .io import load_xenium, materialize_input
from .plotting import make_figures
from .report import build_html_report


def main():
    p = argparse.ArgumentParser(description="Run Xenium QC and spatial application workflow")
    p.add_argument("--config", required=True)
    args = p.parse_args()
    with open(args.config, encoding="utf-8") as f: cfg = yaml.safe_load(f)
    out = Path(cfg["output_dir"]); out.mkdir(parents=True, exist_ok=True)
    root = materialize_input(cfg["input_zip"], out / "input_cache")
    result = analyze(*load_xenium(root), cfg)
    make_figures(result, out, cfg["dataset_kind"])
    report = build_html_report(result, out, cfg)
    print(f"Completed: {out}\nReport: {report}")


if __name__ == "__main__": main()
