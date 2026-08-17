from __future__ import annotations

import base64
import html
import json
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd


def _data_uri(path: Path) -> str:
    return "data:image/png;base64," + base64.b64encode(path.read_bytes()).decode("ascii")


def _table(frame: pd.DataFrame, n: int = 15) -> str:
    return frame.head(n).to_html(index=False, border=0, classes="data", float_format=lambda x: f"{x:.3f}")


def build_html_report(result, outdir, cfg) -> Path:
    out = Path(outdir)
    q = result.cells[result.cells["pass_qc"] == True]
    summary = result.qc_summary
    counts = q["cell_type"].value_counts().rename_axis("cell_type").reset_index(name="cells")
    counts["fraction"] = counts["cells"] / counts["cells"].sum()
    captions = {
        "01_cell_qc.png": "Cell-level yield, detected-feature and control-burden distributions.",
        "02_spatial_qc.png": "Spatial yield reveals localized chemistry, tissue-edge and segmentation effects.",
        "03_cell_types.png": "Conservative marker-module annotation across the tissue section.",
        "04_neighborhood_enrichment.png": "Permutation z-scores quantify preferential spatial adjacency.",
        "05_morans_i.png": "High Moran's I identifies spatially coherent tissue programs.",
        "06_embedding.png": "Expression embedding provides a non-spatial annotation sanity check.",
        "07_myoepithelial_proximity.png": "Distance is a reproducible proxy for epithelial compartment organization.",
        "08_annotation_confidence.png": "Scores and margins expose uncertain labels instead of forcing every cell.",
        "09_cell_type_composition.png": "QC-passing cell-type composition.",
        "10_segmentation_qc.png": "Area-yield coupling and nucleus/cell geometry expose segmentation failure modes.",
    }
    figures = {}
    for path in sorted((out / "figures").glob("*.png")):
        figures[path.name] = f'<figure><img src="{_data_uri(path)}" alt="{html.escape(path.stem)}"><figcaption>{captions.get(path.name, path.stem)}</figcaption></figure>'
    provenance = {
        "dataset": cfg.get("project_name"), "source": cfg.get("dataset_url"),
        "input": cfg.get("input_zip"), "generated_utc": datetime.now(timezone.utc).isoformat(),
        "random_seed": cfg.get("random_seed"), "configuration": cfg,
    }
    cards = [
        ("Cells", f'{summary["cells_total"]:,}'), ("QC pass", f'{summary["pass_rate"]:.1%}'),
        ("Median transcripts", f'{summary["median_transcripts"]:,.0f}'), ("Median genes", f'{summary["median_features"]:,.0f}'),
        ("Panel genes", f'{summary["genes_in_panel"]:,}'), ("Q20 transcripts", f'{summary["q20_transcript_fraction"]:.1%}'),
        ("Assigned (Q20)", f'{summary["q20_assignment_rate"]:.1%}'), ("Unresolved", f'{summary["unresolved_fraction"]:.1%}'),
    ]
    card_html = "".join(f'<div class="card"><span>{html.escape(k)}</span><strong>{v}</strong></div>' for k, v in cards)
    group = lambda names: "".join(figures[n] for n in names if n in figures)
    vendor = pd.DataFrame([summary.get("xoa_metrics", {})]).T.reset_index()
    vendor.columns = ["XOA metric", "value"]
    report = f'''<!doctype html><html><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>{html.escape(cfg.get("project_name", "Xenium report"))}</title><style>
:root{{--ink:#172033;--muted:#667085;--blue:#1769e0;--paper:#fff;--wash:#f2f5f9;--line:#dce3ed}}*{{box-sizing:border-box}}body{{margin:0;background:var(--wash);color:var(--ink);font:15px/1.55 system-ui,Segoe UI,sans-serif}}header{{padding:54px max(5vw,24px);color:white;background:linear-gradient(125deg,#091e42,#1769e0 62%,#16a085)}}header p{{max-width:850px;font-size:18px}}main{{max-width:1280px;margin:auto;padding:28px}}h2{{margin-top:42px;border-bottom:2px solid var(--line);padding-bottom:8px}}.cards{{display:grid;grid-template-columns:repeat(auto-fit,minmax(150px,1fr));gap:12px;margin-top:-54px}}.card,section,figure{{background:white;border:1px solid var(--line);border-radius:12px;box-shadow:0 4px 18px #1720330d}}.card{{padding:16px}}.card span{{display:block;color:var(--muted);font-size:12px;text-transform:uppercase}}.card strong{{font-size:24px}}section{{padding:22px;margin:18px 0}}.grid{{display:grid;grid-template-columns:repeat(auto-fit,minmax(430px,1fr));gap:18px}}figure{{margin:0;padding:14px}}img{{width:100%;height:auto}}figcaption{{color:var(--muted);padding:8px 3px}}table.data{{border-collapse:collapse;width:100%}}.data th,.data td{{padding:7px 10px;border-bottom:1px solid var(--line);text-align:left}}.note{{border-left:4px solid #f59e0b;background:#fffbeb;padding:12px 16px}}code{{white-space:pre-wrap}}footer{{color:var(--muted);padding:32px;text-align:center}}@media(max-width:600px){{.grid{{grid-template-columns:1fr}}}}</style></head>
<body><header><h1>FFPE Human Breast Cancer · Xenium Prime 5K</h1><p>Reproducible QC and spatial tumor-microenvironment analysis separating technical acceptance, segmentation diagnostics, annotation confidence and biological interpretation.</p></header><main>
<div class="cards">{card_html}</div>
<h2>Executive interpretation</h2><section><p><b>Technical:</b> {summary["cells_pass_qc"]:,} of {summary["cells_total"]:,} cells passed configured filters. Q20 assignment was {summary["q20_assignment_rate"]:.1%}; median control fraction was {summary["median_control_fraction"]:.3%}. Spatial QC still requires inspection because acceptable global metrics can hide local artifacts.</p><p><b>Biological:</b> marker modules resolve epithelial, immune, stromal and vascular compartments while {summary["unresolved_fraction"]:.1%} of passing cells remain explicitly unresolved. Neighborhood enrichment and Moran's I are hypothesis-generating and require pathology or independent reference validation.</p><div class="note">Marker-module labels are working annotations, not clinical diagnoses. The spatial graph uses {summary["cells_used_for_spatial_graph"]:,} reproducibly sampled cells for scalability; cell QC and annotation cover all cells.</div></section>
<h2>QC dashboard</h2><div class="grid">{group(["01_cell_qc.png","02_spatial_qc.png","08_annotation_confidence.png","10_segmentation_qc.png"])}</div>
<h2>Cell identity and tissue organization</h2><div class="grid">{group(["03_cell_types.png","04_neighborhood_enrichment.png","06_embedding.png","09_cell_type_composition.png"])}</div>
<h2>Spatial gene programs</h2><section><h3>Top genes by Moran's I</h3>{_table(result.gene_stats, 20)}</section><div class="grid">{group(["05_morans_i.png","07_myoepithelial_proximity.png"])}</div>
<h2>Cell-type composition</h2><section>{_table(counts, 30)}</section>
<h2>Orthogonal XOA metrics</h2><section><p>Vendor metrics are retained separately from recomputed metrics to make disagreements visible.</p>{_table(vendor, 30)}</section>
<h2>Methods and reproducibility</h2><section><ol><li>Stream transcript parquet to calculate QV and assignment metrics.</li><li>Exclude controls; compute per-cell genes, transcripts, control fraction and segmentation geometry flags.</li><li>Library-size normalize, log transform, fit SVD on a reproducible subset, then transform and cluster all passing cells.</li><li>Assign conservative marker-module labels using minimum score and margin thresholds.</li><li>Build a physical k-nearest-neighbor graph, permutation-test cell adjacency and calculate Moran's I.</li></ol><details><summary>Machine-readable provenance</summary><pre><code>{html.escape(json.dumps(provenance, indent=2))}</code></pre></details></section>
</main><footer>Generated by xenium-qc-atlas · self-contained HTML</footer></body></html>'''
    target = out / "xenium_prime5k_qc_analysis_report.html"
    target.write_text(report, encoding="utf-8")
    (out / "run_provenance.json").write_text(json.dumps(provenance, indent=2), encoding="utf-8")
    return target
