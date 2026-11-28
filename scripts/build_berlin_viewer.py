#!/usr/bin/env python3
"""
scripts/build_berlin_viewer.py — Reproducibly generates docs/berlin_triptych_viewer.html
directly from the frozen Phase 2C CSVs in results/ and embeds their exact SHA-256 digests.

Usage:
    uv run python scripts/build_berlin_viewer.py
"""

from __future__ import annotations

import ast
import csv
import hashlib
import json
import pathlib
import numpy as np
from scipy.stats import spearmanr

ROOT = pathlib.Path(__file__).resolve().parent.parent
RESULTS = ROOT / "results"
OUT = ROOT / "docs" / "berlin_triptych_viewer.html"

SOURCE_FILES = [
    "phase2c_berlin_triptych_cells.csv",
    "phase2c_part2_va63038_check.csv",
    "phase2c_reversed_duplicate_kh0121_kh0473.csv",
]


def sha256_file(path: pathlib.Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def build_data() -> dict:
    checksums = {name: sha256_file(RESULTS / name) for name in SOURCE_FILES}

    with (RESULTS / "phase2c_berlin_triptych_cells.csv").open(encoding="utf-8") as f:
        cells = list(csv.DictReader(f))

    with (RESULTS / "phase2c_part2_va63038_check.csv").open(encoding="utf-8") as f:
        part2_raw = list(csv.DictReader(f))

    with (RESULTS / "phase2c_reversed_duplicate_kh0121_kh0473.csv").open(encoding="utf-8") as f:
        revdup_raw = list(csv.DictReader(f))

    def clean_cell(r: dict) -> dict:
        return {
            "half": r["triptych_half"],
            "g": int(r["group_idx"]),
            "s": int(r["slot_idx"]),
            "ap": int(r["master_pnum"]),
            "av": int(r["master_val"]),
            "ac": r["master_col"],
            "al": float(r["master_len_cm"]) if r["master_len_cm"] else None,
            "as": int(r["master_n_subs"]),
            "bp": int(r["copy_pnum"]),
            "bv": int(r["copy_val"]),
            "bc": r["copy_col"],
            "bl": float(r["copy_len_cm"]) if r["copy_len_cm"] else None,
            "bs": int(r["copy_n_subs"]),
            "d": int(r["diff"]),
            "ex": int(r["exact_match"]),
        }

    part1 = [clean_cell(r) for r in cells if r["triptych_half"] == "Part1_vs_UR233"]
    part3 = [clean_cell(r) for r in cells if r["triptych_half"] == "Part3_vs_UR232"]

    # Derive exact unaligned spacer cord ordinals per group directly from bp gaps
    # Part 1 (Khipu B / UR233, p1..p126): p1 unaligned (paired slot G1 S1 absent in Khipu A);
    # spacers at (p-1)%6 == 2 -> p3, p9, p15, ..., p123 (between Slot 2 and Slot 3).
    p1_spacers = {g: 6 * (g - 1) + 3 for g in range(1, 22)}

    # Part 3 (Khipu C / UR232, p1..p125): Groups 1..18 spacers at 6*(g-1)+4 (p4, p10, ..., p106);
    # Group 19 has no spacer; Group 20 spacer at p117; Group 21 spacer at p123 (all between Slot 3 and Slot 4).
    p3_spacers = {g: 6 * (g - 1) + 4 for g in range(1, 19)}
    p3_spacers[20] = 117
    p3_spacers[21] = 123

    # Enrich Part 2 rows with Spearman rho and relative errors computed from exact arrays
    part2 = []
    for r in part2_raw:
        row = dict(r)
        if "Part3_setsum" in r["check"]:
            a = np.array(ast.literal_eval(r["part2_vals"]), dtype=float)
            b = np.array(ast.literal_eval(r["comparison_vals"]), dtype=float)
            row["spearman_rho"] = round(float(spearmanr(a, b).statistic), 3)
            row["sum_rel_err_pct"] = round(float((a.sum() - b.sum()) / b.sum() * 100.0), 1)
            row["mean_abs_rel_err_pct"] = round(float(np.mean(np.abs(a - b) / b) * 100.0), 1)
            row["is_internal"] = True
        else:
            row["spearman_rho"] = None
            row["sum_rel_err_pct"] = None
            row["mean_abs_rel_err_pct"] = None
            row["is_internal"] = False
        part2.append(row)

    revdup = [
        {
            "idx": int(r["aligned_idx"]),
            "o1": int(r["kh0121_ord"]),
            "v1": int(r["kh0121_val"]),
            "l1": float(r["kh0121_len_cm"]),
            "o2": int(r["kh0473_ord_rev"]),
            "v2": int(r["kh0473_val_rev"]),
            "l2": float(r["kh0473_len_cm_rev"]),
            "ex": int(r["val_exact"]),
            "dl": float(r["len_abs_diff_cm"]),
            "unaligned": False,
        }
        for r in revdup_raw
    ]
    # Insert the explicit unaligned cord KH0473 p13 (broken 0-cm stub between p14 and p12)
    revdup_with_unaligned = []
    for item in revdup:
        revdup_with_unaligned.append(item)
        if item["idx"] == 17:
            revdup_with_unaligned.append(
                {
                    "idx": None,
                    "o1": None,
                    "v1": None,
                    "l1": None,
                    "o2": 13,
                    "v2": 0,
                    "l2": 0.0,
                    "ex": 0,
                    "dl": None,
                    "unaligned": True,
                }
            )

    def summarize(rows: list[dict]) -> dict:
        total = len(rows)
        exact_all = sum(c["ex"] for c in rows)
        zero_exact = sum(1 for c in rows if c["ex"] == 1 and c["av"] == 0)
        nz_cells = [c for c in rows if c["av"] > 0 or c["bv"] > 0]
        nz_exact = sum(1 for c in nz_cells if c["ex"] == 1)
        return {
            "total": total,
            "exact_all": exact_all,
            "pct_all": round(100.0 * exact_all / total, 1),
            "zero_exact": zero_exact,
            "nz_total": len(nz_cells),
            "nz_exact": nz_exact,
            "pct_nz": round(100.0 * nz_exact / len(nz_cells), 1),
        }

    return {
        "checksums": checksums,
        "summary": {"part1": summarize(part1), "part3": summarize(part3)},
        "part1": part1,
        "part3": part3,
        "spacers": {"part1": p1_spacers, "part3": p3_spacers},
        "part2": part2,
        "revdup": revdup_with_unaligned,
    }


def render_html(data: dict) -> str:
    s1 = data["summary"]["part1"]
    s3 = data["summary"]["part3"]
    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8" />
<meta name="viewport" content="width=device-width, initial-scale=1.0" />
<title>Putative Berlin Split-Ledger Khipu Triptych — Interactive Cell Viewer</title>
<style>
  :root {{
    --bg: #f6f3ec;
    --surface: #ffffff;
    --surface-alt: #faf8f4;
    --ink: #1e2328;
    --muted: #555e67;
    --border: #dfd9cc;
    --accent: #1d5b54;
    --accent-soft: #e6f2f0;
    --exact-bg: #dff2ec;
    --exact-bd: #6ab5a4;
    --exact-ink: #11473e;
    --near-bg: #fff3da;
    --near-bd: #e0a742;
    --near-ink: #7a4c00;
    --far-bg: #fde4e1;
    --far-bd: #d6685b;
    --far-ink: #7f1d1d;
    --zero-bg: #eef2f3;
    --zero-bd: #c2ccd0;
    --zero-ink: #4d5c66;
    --warn-bg: #fff8e6;
    --warn-bd: #d9a441;
    --shift-ring: #9b2c2c;
    --mono: ui-monospace, SFMono-Regular, Menlo, Consolas, monospace;
    --sans: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Inter, Helvetica, Arial, sans-serif;
    --serif: Georgia, "Times New Roman", serif;
  }}
  * {{ box-sizing: border-box; }}
  body {{
    margin: 0;
    background: var(--bg);
    color: var(--ink);
    font-family: var(--sans);
    line-height: 1.5;
  }}
  .wrap {{
    max-width: 1340px;
    margin: 0 auto;
    padding: 24px 20px 56px;
  }}
  header {{
    background: var(--surface);
    border: 1px solid var(--border);
    border-radius: 14px;
    padding: 24px 28px;
    box-shadow: 0 2px 10px rgba(25, 22, 16, 0.04);
    margin-bottom: 18px;
  }}
  .eyebrow {{
    display: flex;
    flex-wrap: wrap;
    gap: 8px;
    align-items: center;
    margin-bottom: 10px;
    font-size: 12px;
    font-weight: 600;
    letter-spacing: 0.04em;
    text-transform: uppercase;
  }}
  .pill {{
    display: inline-flex;
    align-items: center;
    padding: 3px 10px;
    border-radius: 999px;
    background: var(--accent-soft);
    color: var(--accent);
    border: 1px solid #b8dcd6;
  }}
  .pill.warn {{
    background: var(--warn-bg);
    color: #784f00;
    border-color: var(--warn-bd);
  }}
  .pill.alert {{
    background: var(--far-bg);
    color: var(--far-ink);
    border-color: var(--far-bd);
  }}
  h1 {{
    margin: 0 0 8px;
    font-family: var(--serif);
    font-size: clamp(22px, 3vw, 31px);
    line-height: 1.2;
    color: #161b20;
  }}
  .subtitle {{
    margin: 0 0 14px;
    color: var(--muted);
    font-size: 14.5px;
    max-width: 1020px;
  }}
  .meta-links {{
    display: flex;
    flex-wrap: wrap;
    gap: 14px;
    font-size: 13px;
    color: var(--muted);
  }}
  .meta-links a {{
    color: var(--accent);
    text-decoration: none;
    font-weight: 600;
  }}
  .meta-links a:hover {{ text-decoration: underline; }}

  /* Non-directional 3-piece summary cards */
  .arch {{
    display: grid;
    grid-template-columns: 1.2fr 1fr 1fr;
    gap: 12px;
    margin-bottom: 12px;
  }}
  @media (max-width: 940px) {{
    .arch {{ grid-template-columns: 1fr; }}
  }}
  .card {{
    background: var(--surface);
    border: 1px solid var(--border);
    border-radius: 12px;
    padding: 16px 18px;
    box-shadow: 0 1px 6px rgba(25, 22, 16, 0.03);
  }}
  .card h3 {{
    margin: 0 0 6px;
    font-size: 14px;
    display: flex;
    justify-content: space-between;
    align-items: center;
    gap: 8px;
    flex-wrap: wrap;
  }}
  .tag {{
    font-family: var(--mono);
    font-size: 11.5px;
    padding: 2px 7px;
    border-radius: 6px;
    background: var(--surface-alt);
    border: 1px solid var(--border);
    color: var(--muted);
  }}
  .card p {{
    margin: 4px 0;
    font-size: 12.8px;
    color: var(--muted);
  }}
  .stat-row {{
    display: flex;
    flex-wrap: wrap;
    gap: 7px;
    margin-top: 8px;
    font-family: var(--mono);
    font-size: 11.8px;
  }}
  .stat-chip {{
    background: var(--surface-alt);
    border: 1px solid var(--border);
    padding: 3px 8px;
    border-radius: 6px;
  }}
  .stat-chip.strong {{
    background: var(--exact-bg);
    border-color: var(--exact-bd);
    color: var(--exact-ink);
    font-weight: 700;
  }}
  .neutral-banner {{
    background: var(--surface);
    border-left: 4px solid var(--accent);
    border-top: 1px solid var(--border);
    border-right: 1px solid var(--border);
    border-bottom: 1px solid var(--border);
    border-radius: 8px;
    padding: 10px 14px;
    font-size: 12.8px;
    color: var(--muted);
    margin-bottom: 16px;
  }}

  /* Tabs */
  .tabs {{
    display: flex;
    flex-wrap: wrap;
    gap: 8px;
    margin-bottom: 14px;
  }}
  .tab-btn {{
    appearance: none;
    border: 1px solid var(--border);
    background: var(--surface);
    color: var(--ink);
    padding: 10px 14px;
    border-radius: 10px;
    font-size: 13px;
    font-weight: 600;
    cursor: pointer;
    transition: all 0.15s ease;
  }}
  .tab-btn:hover {{ border-color: var(--accent); }}
  .tab-btn.active {{
    background: var(--accent);
    color: #ffffff;
    border-color: var(--accent);
  }}

  /* Workspace */
  .workspace {{
    display: grid;
    grid-template-columns: minmax(0, 1fr) 340px;
    gap: 16px;
    align-items: start;
  }}
  @media (max-width: 1060px) {{
    .workspace {{ grid-template-columns: 1fr; }}
  }}
  .toolbar {{
    display: flex;
    flex-wrap: wrap;
    justify-content: space-between;
    gap: 10px;
    align-items: center;
    margin-bottom: 10px;
    padding-bottom: 10px;
    border-bottom: 1px solid var(--border);
  }}
  .btn-group {{
    display: flex;
    flex-wrap: wrap;
    gap: 5px;
    align-items: center;
  }}
  .btn-group label {{
    font-size: 11.5px;
    font-weight: 700;
    color: var(--muted);
    margin-right: 3px;
    text-transform: uppercase;
    letter-spacing: 0.04em;
  }}
  .ctrl {{
    border: 1px solid var(--border);
    background: var(--surface-alt);
    color: var(--ink);
    padding: 5px 9px;
    border-radius: 7px;
    font-size: 12px;
    font-weight: 600;
    cursor: pointer;
  }}
  .ctrl.active {{
    background: #1e293b;
    color: #fff;
    border-color: #1e293b;
  }}
  .ctrl.highlight-btn {{
    background: #fff3da;
    border-color: #d9a441;
    color: #6b4200;
  }}
  .ctrl.highlight-btn.active {{
    background: #9b2c2c;
    border-color: #9b2c2c;
    color: #fff;
  }}

  .callout {{
    border-radius: 9px;
    padding: 10px 13px;
    font-size: 12.6px;
    margin-bottom: 12px;
    border: 1px solid var(--border);
    background: var(--surface-alt);
    color: var(--ink);
  }}
  .callout.pattern {{
    background: #fff9eb;
    border-color: #e2b857;
    color: #5c3c00;
  }}
  .callout.color-note {{
    background: #f1f6f9;
    border-color: #b7d0df;
    color: #1c3d52;
  }}

  .legend {{
    display: flex;
    flex-wrap: wrap;
    gap: 12px;
    font-size: 11.8px;
    color: var(--muted);
    margin-bottom: 10px;
  }}
  .leg-item {{
    display: inline-flex;
    align-items: center;
    gap: 5px;
  }}
  .swatch {{
    width: 13px;
    height: 13px;
    border-radius: 3px;
    border: 1px solid transparent;
    display: inline-block;
  }}
  .swatch.unaligned-sw {{
    background: repeating-linear-gradient(45deg, #f0ece2, #f0ece2 4px, #e2dccd 4px, #e2dccd 8px);
    border-color: #bdb5a3;
  }}

  /* Matrix */
  .table-scroll {{ overflow-x: auto; }}
  table.matrix {{
    width: 100%;
    border-collapse: separate;
    border-spacing: 5px;
    font-family: var(--mono);
    font-size: 12px;
  }}
  table.matrix th {{
    font-family: var(--sans);
    font-size: 11.8px;
    color: var(--muted);
    font-weight: 700;
    text-align: center;
    padding: 5px 4px;
  }}
  table.matrix th.row-hdr {{
    text-align: right;
    padding-right: 7px;
    width: 72px;
    white-space: nowrap;
  }}
  table.matrix th.spacer-hdr {{
    color: #7c7464;
    font-style: italic;
    font-weight: 600;
    font-size: 11px;
  }}
  td.cell {{
    border-radius: 8px;
    padding: 6px 5px;
    text-align: center;
    cursor: pointer;
    border: 1.5px solid transparent;
    transition: transform 0.1s ease, box-shadow 0.12s ease, opacity 0.15s ease;
    min-width: 96px;
    vertical-align: middle;
  }}
  td.cell:hover, td.cell:focus-visible {{
    transform: translateY(-1px);
    box-shadow: 0 3px 8px rgba(0, 0, 0, 0.12);
    outline: 2px solid var(--accent);
    outline-offset: 1px;
  }}
  td.cell.selected {{
    outline: 2.5px solid #161b20;
    outline-offset: 1px;
  }}
  td.cell.pattern-ring {{
    box-shadow: inset 0 0 0 2px var(--shift-ring);
  }}
  td.cell.dim {{ opacity: 0.2; }}
  td.cell.exact {{
    background: var(--exact-bg);
    border-color: var(--exact-bd);
    color: var(--exact-ink);
  }}
  td.cell.exact-zero {{
    background: var(--zero-bg);
    border-color: var(--zero-bd);
    color: var(--zero-ink);
  }}
  td.cell.near {{
    background: var(--near-bg);
    border-color: var(--near-bd);
    color: var(--near-ink);
  }}
  td.cell.far {{
    background: var(--far-bg);
    border-color: var(--far-bd);
    color: var(--far-ink);
  }}
  td.cell.unaligned {{
    background: repeating-linear-gradient(45deg, #f3f0e8, #f3f0e8 6px, #e9e3d5 6px, #e9e3d5 12px);
    border-color: #cbc3b2;
    color: #6b6558;
    min-width: 78px;
    cursor: default;
  }}
  td.cell.unaligned:hover {{
    transform: none;
    box-shadow: none;
  }}
  .c-main {{
    font-weight: 700;
    font-size: 12.5px;
    display: block;
  }}
  .c-sub {{
    font-size: 10.3px;
    opacity: 0.86;
    display: block;
    margin-top: 2px;
  }}
  .fiber-dot {{
    display: inline-block;
    width: 9px;
    height: 9px;
    border-radius: 50%;
    border: 1px solid rgba(0,0,0,0.28);
    vertical-align: middle;
    margin-right: 3px;
  }}
  tfoot tr td {{
    text-align: center;
    font-family: var(--sans);
    font-size: 11.2px;
    color: var(--muted);
    padding: 3px 2px;
    font-weight: 600;
  }}
  tfoot tr.nz-row td {{
    color: var(--exact-ink);
    font-weight: 700;
    padding-top: 8px;
  }}

  /* Inspector */
  .inspector {{
    position: sticky;
    top: 16px;
  }}
  .inspector h3 {{
    margin-bottom: 8px;
    font-size: 15px;
  }}
  .kv-grid {{
    display: grid;
    grid-template-columns: 1fr 1fr;
    gap: 8px;
    margin: 10px 0;
  }}
  .kv-box {{
    background: var(--surface-alt);
    border: 1px solid var(--border);
    border-radius: 9px;
    padding: 9px 10px;
  }}
  .kv-box .k {{
    font-size: 10.8px;
    color: var(--muted);
    text-transform: uppercase;
    letter-spacing: 0.04em;
    font-weight: 700;
  }}
  .kv-box .v {{
    font-family: var(--mono);
    font-size: 16px;
    font-weight: 700;
    margin-top: 2px;
  }}
  .kv-box .s {{
    font-size: 11.3px;
    color: var(--muted);
    margin-top: 2px;
  }}
  .bar-track {{
    background: #ebe6db;
    height: 9px;
    border-radius: 999px;
    overflow: hidden;
    margin: 4px 0 8px;
  }}
  .bar-fill {{
    height: 100%;
    border-radius: 999px;
  }}
  .note-box {{
    margin-top: 12px;
    background: var(--warn-bg);
    border: 1px solid var(--warn-bd);
    border-radius: 9px;
    padding: 10px 12px;
    font-size: 12px;
    color: #5e3e00;
  }}

  /* Data tables */
  table.clean {{
    width: 100%;
    border-collapse: collapse;
    font-size: 12.8px;
  }}
  table.clean th, table.clean td {{
    padding: 8px 10px;
    border-bottom: 1px solid var(--border);
    text-align: left;
  }}
  table.clean th {{
    background: var(--surface-alt);
    font-size: 11.5px;
    text-transform: uppercase;
    letter-spacing: 0.03em;
    color: var(--muted);
  }}
  table.clean td.num {{
    font-family: var(--mono);
    text-align: right;
  }}
  tr.museum-flag {{
    background: #fff4ec;
  }}
  tr.unaligned-row {{
    background: #f3f0e8;
    color: #5d5649;
  }}
  .mini-bars {{
    display: flex;
    flex-direction: column;
    gap: 3px;
    min-width: 140px;
  }}
  .mini-bar-row {{
    display: flex;
    align-items: center;
    gap: 6px;
    font-family: var(--mono);
    font-size: 11px;
  }}
  .mini-bar {{
    height: 7px;
    border-radius: 4px;
  }}
  .checklist {{
    display: grid;
    grid-template-columns: repeat(3, 1fr);
    gap: 10px;
    margin: 12px 0 16px;
  }}
  @media (max-width: 860px) {{
    .checklist {{ grid-template-columns: 1fr; }}
  }}
  .chk-card {{
    background: #fff6f2;
    border: 1px solid #e4a38f;
    border-radius: 9px;
    padding: 10px 12px;
    font-size: 12.4px;
  }}
  .chk-card strong {{
    color: #7f1d1d;
    display: block;
    margin-bottom: 3px;
  }}

  details.hashes {{
    margin-top: 18px;
    background: var(--surface);
    border: 1px solid var(--border);
    border-radius: 10px;
    padding: 10px 14px;
    font-size: 12.3px;
    color: var(--muted);
  }}
  details.hashes summary {{
    cursor: pointer;
    font-weight: 700;
    color: var(--ink);
  }}
  details.hashes code {{
    font-family: var(--mono);
    font-size: 11.3px;
    word-break: break-all;
  }}
  .hidden {{ display: none !important; }}
  footer {{
    margin-top: 24px;
    font-size: 12.3px;
    color: var(--muted);
    text-align: center;
    line-height: 1.7;
  }}
  footer a {{ color: var(--accent); text-decoration: none; font-weight: 600; }}
</style>
</head>
<body>
<div class="wrap">
  <header>
    <div class="eyebrow">
      <span class="pill">Phase 2C · Exploratory Analysis</span>
      <span class="pill warn">Putative Links — Pending Museum Verification</span>
    </div>
    <h1>Putative Berlin Split-Ledger Khipu Triptych &amp; Reversed Physical Duplicate</h1>
    <p class="subtitle">
      Interactive cell-by-cell explorer for the descriptive <code>611 × 611</code> rare-number sweep of OKR v2.1.0
      (<code>538</code> pre-filtered pairs from ~186k possible pairs). Inspect aligned cells, unaligned structural spacers,
      non-zero match rates, fiber color-code contrasts, and physical cord-length forensics across <strong>Khipu A</strong>
      (<code>VA 42518</code>), <strong>Khipu B</strong> (<code>VA 63038 b</code>), <strong>Khipu C</strong> (<code>VA 63038(B)</code>),
      and the putative reversed duplicate <strong>KH0121 ↔ KH0473</strong>.
    </p>
    <div class="meta-links">
      <span><strong>Author:</strong> Maria Contreras (Disentangle Lab)</span>
      <a href="https://orcid.org/0009-0003-0749-4971" target="_blank" rel="noopener">ORCID 0009-0003-0749-4971</a>
      <a href="https://doi.org/10.5281/zenodo.24610674" target="_blank" rel="noopener">Zenodo Concept DOI 10.5281/zenodo.24610674</a>
      <a href="https://github.com/Disentangle-Lab/khipu-semantic-alignment/releases/tag/v1.0.1" target="_blank" rel="noopener">GitHub Release v1.0.1</a>
    </div>
  </header>

  <section class="arch" aria-label="Three-piece putative ledger overview">
    <div class="card">
      <h3>Khipu A · Longer 45-Group Ledger <span class="tag">KH0192 / VA 42518</span></h3>
      <p>Cataloged as <strong>Pachacámac</strong> (<code>AS175 / UR1175</code>). 45 five-cord groups: <strong>Part 1</strong> (G1–21, 104 cords; G1 S1 absent), <strong>Part 2</strong> (G22–24, 15 cords), <strong>Part 3</strong> (G25–45, 105 cords).</p>
      <div class="stat-row">
        <span class="stat-chip">k = 5</span>
        <span class="stat-chip">Med len: 65.0 cm</span>
        <span class="stat-chip">114 subsidiaries</span>
      </div>
    </div>
    <div class="card">
      <h3>Khipu B · Shorter Part 1 Match <span class="tag">KH0470 / VA 63038 b</span></h3>
      <p>Cataloged as <strong>Huacho</strong> (<code>UR233</code>). 126 pendants = 104 aligned + <code>p1</code> unaligned + 21 zero spacers. <strong>Sweep Rank #1 / 538</strong> (36 exact rare matches vs. max background 8).</p>
      <div class="stat-row">
        <span class="stat-chip strong">Non-zero: {s1['nz_exact']}/{s1['nz_total']} ({s1['pct_nz']}%)</span>
        <span class="stat-chip">All: {s1['exact_all']}/{s1['total']} ({s1['pct_all']}%, incl. {s1['zero_exact']} 0=0)</span>
        <span class="stat-chip">r = 0.9823 · Med 40.0 cm</span>
      </div>
    </div>
    <div class="card">
      <h3>Khipu C · Shorter Part 3 Match <span class="tag">KH0469 / VA 63038(B)</span></h3>
      <p>Cataloged as <strong>Huacho</strong> (<code>UR232</code>). 125 pendants = 105 aligned + 20 zero spacers. Surfaced via targeted inspection of Berlin lot <code>VA 63038</code> (outside ±10 sweep lag).</p>
      <div class="stat-row">
        <span class="stat-chip strong">Non-zero: {s3['nz_exact']}/{s3['nz_total']} ({s3['pct_nz']}%)</span>
        <span class="stat-chip">All: {s3['exact_all']}/{s3['total']} ({s3['pct_all']}%, incl. {s3['zero_exact']} 0=0)</span>
        <span class="stat-chip">r = 0.9446 · Med 51.5 cm</span>
      </div>
    </div>
  </section>

  <div class="neutral-banner">
    <strong>Non-directional labeling:</strong> Labels <strong>Khipu A</strong> (longer), <strong>Khipu B</strong>, and <strong>Khipu C</strong> (shorter) avoid assuming copying direction.
    Numerical correspondence alone cannot determine whether B and C were copied from A, A compiled from B and C, or all three transcribed from a shared third source, pending Ethnologisches Museum Berlin verification.
  </div>

  <nav class="tabs" role="tablist">
    <button class="tab-btn active" data-tab="part1">Part 1: Khipu A (G1–21) ↔ Khipu B ({s1['nz_exact']}/{s1['nz_total']} non-zero · {s1['exact_all']}/{s1['total']} all)</button>
    <button class="tab-btn" data-tab="part3">Part 3: Khipu A (G25–45) ↔ Khipu C ({s3['nz_exact']}/{s3['nz_total']} non-zero · {s3['exact_all']}/{s3['total']} all)</button>
    <button class="tab-btn" data-tab="part2">Part 2: Internal Bridge (G22–24) &amp; Lot Absence</button>
    <button class="tab-btn" data-tab="revdup">Reversed Duplicate: VA 24370(B) ↔ VA 66834 (27/29)</button>
  </nav>

  <!-- MATRIX VIEW (Part 1 & Part 3) -->
  <section id="matrixView" class="workspace">
    <div class="card">
      <div class="toolbar">
        <div class="btn-group" id="filterBtns">
          <label>Filter:</label>
          <button class="ctrl active" data-filter="all">All cells</button>
          <button class="ctrl" data-filter="nz-exact">Non-zero exact</button>
          <button class="ctrl" data-filter="exact">All exact (incl. 0=0)</button>
          <button class="ctrl" data-filter="diff">Adjusted (Δ ≠ 0)</button>
        </div>
        <div class="btn-group" id="modeBtns">
          <label>Show:</label>
          <button class="ctrl active" data-mode="values">Values</button>
          <button class="ctrl" data-mode="diff">Δ Diff</button>
          <button class="ctrl" data-mode="lengths">Lengths (cm)</button>
          <button class="ctrl" data-mode="colors">Colors</button>
          <button class="ctrl" data-mode="subs">Subsidiaries</button>
        </div>
      </div>

      <div class="toolbar" style="border-bottom:none;margin-bottom:6px;padding-bottom:0;">
        <div class="btn-group">
          <button class="ctrl active" id="toggleSpacers" aria-pressed="true">Show unaligned cords &amp; spacers</button>
          <button class="ctrl highlight-btn" id="toggleSlot4Pattern" aria-pressed="false">Highlight Slot 4 G1–G6 systematic shift (+10 to +20)</button>
        </div>
      </div>

      <div id="patternCallout" class="callout pattern">
        <strong>Structured discrepancy in Part 1, Slot 4 (Groups 1–6):</strong> Rather than random noise, Khipu B exceeds Khipu A by a systematic <code>+10</code> to <code>+20</code> increment across the first six groups
        (<code>100→110</code> in G1, G2, G4, G5; <code>100→120</code> in G3; <code>58→70</code> in G6), then switches to exact parity in G7–G19 (except <code>10→100</code> in G20). Click <em>Highlight Slot 4 G1–G6</em> to isolate this block.
      </div>

      <div id="colorCallout" class="callout color-note hidden">
        <strong>Color-code shift caveat (cataloger vocabulary vs. physical dye):</strong> Per-slot color codes shift systematically between Khipu A and Khipus B/C (e.g., <code>LB→RB/PR</code>, <code>YG→GG/BG</code>, <code>YB→AB</code>).
        Because Khipu A was originally documented under Ascher conventions (<code>AS175</code>) whereas B and C were recorded under Urton/KFG coding (<code>UR233</code>, <code>UR232</code>), these shifts may reflect <em>inter-cataloger color taxonomy</em>, genuine physical dye differences, or both—resolvable only by side-by-side museum inspection.
      </div>

      <div class="legend">
        <span class="leg-item"><span class="swatch" style="background:var(--exact-bg);border-color:var(--exact-bd);"></span>Non-zero exact</span>
        <span class="leg-item"><span class="swatch" style="background:var(--zero-bg);border-color:var(--zero-bd);"></span>Zero exact (0 = 0)</span>
        <span class="leg-item"><span class="swatch" style="background:var(--near-bg);border-color:var(--near-bd);"></span>Small shift (1 ≤ |Δ| ≤ 10)</span>
        <span class="leg-item"><span class="swatch" style="background:var(--far-bg);border-color:var(--far-bd);"></span>Larger shift (|Δ| &gt; 10)</span>
        <span class="leg-item"><span class="swatch unaligned-sw"></span>Unaligned cord / 0-spacer (excluded from denominator)</span>
      </div>

      <div class="table-scroll">
        <table class="matrix" id="matrixTable"></table>
      </div>
    </div>

    <aside class="card inspector" id="inspectorPanel">
      <h3>Cell Forensics Inspector</h3>
      <p id="inspHint" style="color:var(--muted);font-size:12.8px;">Hover, click, or press Enter/Space on any aligned cell to inspect cord ordinals, values, fiber colors, physical cord lengths, and subsidiary counts.</p>
      <div id="inspBody" class="hidden"></div>
      <div class="note-box">
        <strong>Methodological note:</strong> Cell-shuffling <code>p ≤ 1e-4</code> values (<code>N_perm = 10,000</code>) are <em>conditional post-selection permutations</em> uncorrected for the ~186k-pair corpus sweep. Denominators reflect aligned data cords only; unaligned cords and structural 0-spacers are shown explicitly in the hatched cells.
      </div>
    </aside>
  </section>

  <!-- PART 2 VIEW -->
  <section id="part2View" class="card hidden">
    <h3 style="font-size:17px;margin-bottom:8px;">Part 2 Verification: Internal Summary Bridge (Groups 22–24) &amp; Absence from Lot VA 63038</h3>
    <p style="font-size:13.3px;color:var(--muted);margin-bottom:14px;">
      Between Part 1 (Groups 1–21) and Part 3 (Groups 25–45), Khipu A (<code>VA 42518</code>) contains three five-cord groups (Groups 22–24, 15 cords).
      Internal comparison is <strong>consistent with Ascher's Note 5</strong>: all three 5-slot rows show perfect rank-order agreement (Spearman <code>ρ = 1.000</code>, Pearson <code>r = 0.978–0.995</code>) with Part 3's column-by-column 7-group set sums,
      though <em>exact cell matches are only 2/5, 2/5, and 1/5</em> and row sums differ by <code>−2.9%</code> to <code>−16.0%</code> (high Pearson <code>r</code> with <code>n = 5</code> is strongly leveraged by large values 200–445).
      Sliding-window checks across all three khipus in Berlin lot <strong>VA 63038</strong> (<code>UR231</code>, <code>UR232</code>, <code>UR233</code>) find no Part 2 counterpart (max <code>2/15</code> exact matches).
    </p>
    <div class="table-scroll">
      <table class="clean" id="part2Table"></table>
    </div>
  </section>

  <!-- REVERSED DUPLICATE VIEW -->
  <section id="revdupView" class="card hidden">
    <h3 style="font-size:17px;margin-bottom:8px;">Putative Reversed Physical Duplicate: KH0121 (VA 24370(B), Márquez) ↔ KH0473 (VA 66834, unknown)</h3>
    <p style="font-size:13.3px;color:var(--muted);margin-bottom:10px;">
      Ranked <strong>#3 / 538</strong> in the descriptive <code>611 × 611</code> sweep (<code>ord_rare_ex = 9</code>, <code>rare_dice = 0.900</code>). Reading <strong>KH0473 / UR237</strong> in reverse
      matches <strong>KH0121 / UR1107</strong> on <strong>27 / 29 (93.1%)</strong> aligned pendant values (<code>ρ = 0.9332</code>) and near-identical cord lengths
      (<code>r = 0.9767</code>, median <code>|Δlen| = 1.5 cm</code>, <code>66.5</code> vs. <code>66.2 cm</code>).
    </p>
    <div class="callout">
      <strong>Competing hypotheses for Ethnologisches Museum Berlin:</strong>
      (1) <em>Same physical artifact</em> cataloged twice in opposite directions (which would restore Márquez provenance to <code>VA 66834</code>), <strong>or</strong>
      (2) <em>Two distinct companion copies</em> intentionally constructed to the same standard cord length (~66 cm). The three highlighted anomalies below are the highest-leverage physical checkpoints to resolve (1) vs. (2):
    </div>
    <div class="checklist">
      <div class="chk-card">
        <strong>Checkpoint 1 · Aligned #5 (28 ↔ 820)</strong>
        <code>KH0121 p5 = 28</code> (57.5 cm) vs. <code>KH0473 p26 = 820</code> (57.0 cm, <code>|Δlen| = 0.5 cm</code>). Nearly identical cord length with transposed/shifted decade reading—check knot tiers physically.
      </div>
      <div class="chk-card">
        <strong>Checkpoint 2 · Aligned #7 (610 ↔ 66)</strong>
        <code>KH0121 p7 = 610</code> (63.0 cm) vs. <code>KH0473 p24 = 66</code> (57.0 cm, <code>|Δlen| = 6.0 cm</code>). Both knot value and length diverge—check whether a lower knot cluster is missing or re-read.
      </div>
      <div class="chk-card">
        <strong>Checkpoint 3 · Unpaired Cord KH0473 p13</strong>
        <code>KH0473</code> has <strong>30</strong> pendants vs. <strong>29</strong> in <code>KH0121</code>: reversed cord <code>p13</code> (shown below between aligned #17 and #18) is an unaligned <code>0 cm</code> broken stub with no counterpart in <code>KH0121</code>.
      </div>
    </div>
    <div class="table-scroll">
      <table class="clean" id="revdupTable"></table>
    </div>
  </section>

  <details class="hashes">
    <summary>Source CSV SHA-256 Checksums &amp; Reproducible Build Provenance</summary>
    <p style="margin:8px 0 6px;">Generated deterministically by <code>scripts/build_berlin_viewer.py</code> from frozen <code>v1.0.1</code> tables in <code>results/</code> (derived from OKR v2.1.0, DOI <a href="https://doi.org/10.5281/zenodo.18025748" target="_blank" rel="noopener">10.5281/zenodo.18025748</a>):</p>
    <ul style="margin:4px 0 0;padding-left:18px;">
      <li><code>results/phase2c_berlin_triptych_cells.csv</code>: <code>{data['checksums']['phase2c_berlin_triptych_cells.csv']}</code></li>
      <li><code>results/phase2c_part2_va63038_check.csv</code>: <code>{data['checksums']['phase2c_part2_va63038_check.csv']}</code></li>
      <li><code>results/phase2c_reversed_duplicate_kh0121_kh0473.csv</code>: <code>{data['checksums']['phase2c_reversed_duplicate_kh0121_kh0473.csv']}</code></li>
    </ul>
  </details>

  <footer>
    Data derived from the Open Khipu Repository (OKR v2.1.0, <a href="https://doi.org/10.5281/zenodo.18025748" target="_blank" rel="noopener">DOI 10.5281/zenodo.18025748</a>)
    and Khipu Field Guide Excel workbooks (pinned commit <code>167b8df</code>) ·
    Code MIT, Derived Visuals CC-BY-4.0 ·
    <a href="https://github.com/Disentangle-Lab/khipu-semantic-alignment/releases/tag/v1.0.1" target="_blank" rel="noopener">Disentangle-Lab/khipu-semantic-alignment (v1.0.1)</a>
  </footer>
</div>

<script>
const DATA = {json.dumps(data, separators=(",", ":"))};

const FIBER_COLOR = {{
  "W": "#f5f2ea",
  "LB": "#c89f72",
  "RB": "#9c563a",
  "PR": "#8d4739",
  "RL": "#b56952",
  "YG": "#b5a659",
  "GG": "#84915f",
  "BG": "#6e7a5e",
  "YB": "#b58b4c",
  "AB": "#9c7447",
  "PB": "#7d5b3b",
  "LG": "#b9ba9f",
  "MB": "#6a4a33",
  "KB": "#4a3324",
  "HB": "#5a3c29",
  "0G": "#777c72"
}};

function primarySwatch(code) {{
  if (!code) return "#cccccc";
  const base = code.split(/[:(/]/)[0].trim();
  return FIBER_COLOR[base] || "#a89882";
}}

let activeTab = "part1";
let activeFilter = "all";
let activeMode = "values";
let showSpacers = true;
let highlightSlot4 = false;
let selectedCellKey = null;
let lastInspectedKey = null;

function cellClass(c) {{
  if (c.ex === 1) return c.av === 0 ? "exact-zero" : "exact";
  return Math.abs(c.d) <= 10 ? "near" : "far";
}}

function passesFilter(c) {{
  if (highlightSlot4 && activeTab === "part1") {{
    return c.s === 4 && c.g >= 1 && c.g <= 6;
  }}
  if (activeFilter === "all") return true;
  if (activeFilter === "exact") return c.ex === 1;
  if (activeFilter === "nz-exact") return c.ex === 1 && c.av > 0;
  if (activeFilter === "diff") return c.ex === 0;
  return true;
}}

function cellHtml(c) {{
  if (activeMode === "values") {{
    if (c.ex === 1) return `<span class="c-main">${{c.av}}</span><span class="c-sub">exact (= ${{c.bv}})</span>`;
    const sign = c.d > 0 ? `+${{c.d}}` : `${{c.d}}`;
    return `<span class="c-main">${{c.av}} ↔ ${{c.bv}}</span><span class="c-sub">Δ ${{sign}}</span>`;
  }}
  if (activeMode === "diff") {{
    if (c.ex === 1) return `<span class="c-main">0</span><span class="c-sub">val=${{c.av}}</span>`;
    const sign = c.d > 0 ? `+${{c.d}}` : `${{c.d}}`;
    return `<span class="c-main">${{sign}}</span><span class="c-sub">${{c.av}} ↔ ${{c.bv}}</span>`;
  }}
  if (activeMode === "lengths") {{
    const al = c.al !== null ? c.al.toFixed(1) : "—";
    const bl = c.bl !== null ? c.bl.toFixed(1) : "—";
    const dl = (c.al !== null && c.bl !== null) ? (c.bl - c.al).toFixed(1) : "—";
    return `<span class="c-main">${{al}} / ${{bl}}</span><span class="c-sub">cm (Δ ${{dl}})</span>`;
  }}
  if (activeMode === "colors") {{
    const acShort = c.ac.length > 8 ? c.ac.slice(0, 8) + "…" : c.ac;
    const bcShort = c.bc.length > 8 ? c.bc.slice(0, 8) + "…" : c.bc;
    return `<span class="c-main"><span class="fiber-dot" style="background:${{primarySwatch(c.ac)}}"></span>${{acShort}}</span><span class="c-sub"><span class="fiber-dot" style="background:${{primarySwatch(c.bc)}}"></span>${{bcShort}}</span>`;
  }}
  if (activeMode === "subs") {{
    return `<span class="c-main">${{c.as}} / ${{c.bs}}</span><span class="c-sub">A / ${{activeTab === "part1" ? "B" : "C"}} subs</span>`;
  }}
  return "";
}}

function renderInspector(c) {{
  const otherLabel = activeTab === "part1" ? "Khipu B (VA 63038 b)" : "Khipu C (VA 63038(B))";
  const khipuAGroup = activeTab === "part1" ? c.g : c.g + 24;
  const badge = c.ex === 1
    ? (c.av === 0
        ? `<span class="pill" style="background:var(--zero-bg);color:var(--zero-ink);border-color:var(--zero-bd);">Zero Match (0 = 0)</span>`
        : `<span class="pill" style="background:var(--exact-bg);color:var(--exact-ink);border-color:var(--exact-bd);">Non-Zero Exact Match</span>`)
    : `<span class="pill warn">Value Difference (B/C − A = ${{c.d > 0 ? "+" + c.d : c.d}})</span>`;

  const maxLen = Math.max(c.al || 0, c.bl || 0, 85);
  const alPct = c.al ? Math.min(100, (c.al / maxLen) * 100).toFixed(1) : 0;
  const blPct = c.bl ? Math.min(100, (c.bl / maxLen) * 100).toFixed(1) : 0;

  document.getElementById("inspHint").classList.add("hidden");
  const body = document.getElementById("inspBody");
  body.classList.remove("hidden");
  body.innerHTML = `
    <div style="margin-bottom:8px;">${{badge}}</div>
    <div style="font-size:13px;font-weight:700;margin-bottom:6px;">
      Khipu A Group ${{khipuAGroup}} · ${{activeTab === "part1" ? "Khipu B" : "Khipu C"}} Group ${{c.g}} · Slot ${{c.s}}
    </div>
    <div class="kv-grid">
      <div class="kv-box">
        <div class="k">Khipu A (VA 42518)</div>
        <div class="v">${{c.av}}</div>
        <div class="s">Cord <code>p${{c.ap}}</code> · ${{c.as}} sub(s)</div>
        <div class="s"><span class="fiber-dot" style="background:${{primarySwatch(c.ac)}}"></span><code>${{c.ac}}</code></div>
      </div>
      <div class="kv-box">
        <div class="k">${{otherLabel}}</div>
        <div class="v">${{c.bv}}</div>
        <div class="s">Cord <code>p${{c.bp}}</code> · ${{c.bs}} sub(s)</div>
        <div class="s"><span class="fiber-dot" style="background:${{primarySwatch(c.bc)}}"></span><code>${{c.bc}}</code></div>
      </div>
    </div>
    <div style="margin-top:10px;font-size:12px;">
      <strong>Physical Cord Length Comparison</strong>
      <div style="display:flex;justify-content:space-between;margin-top:4px;">
        <span>Khipu A (<code>VA 42518</code>)</span><strong>${{c.al !== null ? c.al.toFixed(1) + " cm" : "n/a"}}</strong>
      </div>
      <div class="bar-track"><div class="bar-fill" style="width:${{alPct}}%;background:#1d5b54;"></div></div>
      <div style="display:flex;justify-content:space-between;">
        <span>${{otherLabel.split(" ")[0] + " " + otherLabel.split(" ")[1]}}</span><strong>${{c.bl !== null ? c.bl.toFixed(1) + " cm" : "n/a"}}</strong>
      </div>
      <div class="bar-track"><div class="bar-fill" style="width:${{blPct}}%;background:#d98a2b;"></div></div>
    </div>
  `;
}}

function renderMatrix() {{
  const rows = DATA[activeTab];
  const byKey = new Map();
  rows.forEach(c => byKey.set(`${{c.g}}-${{c.s}}`, c));
  const spacers = DATA.spacers[activeTab];
  const spacerAfterSlot = activeTab === "part1" ? 2 : 3;

  // Callout visibility
  document.getElementById("patternCallout").classList.toggle("hidden", activeTab !== "part1");
  document.getElementById("toggleSlot4Pattern").classList.toggle("hidden", activeTab !== "part1");
  document.getElementById("colorCallout").classList.toggle("hidden", activeMode !== "colors");

  const slotNzTotal = [0,0,0,0,0];
  const slotNzExact = [0,0,0,0,0];
  const slotZeroExact = [0,0,0,0,0];
  const slotTotal = [0,0,0,0,0];
  const slotAllExact = [0,0,0,0,0];
  rows.forEach(c => {{
    const i = c.s - 1;
    slotTotal[i]++;
    if (c.ex === 1) slotAllExact[i]++;
    if (c.av > 0 || c.bv > 0) {{
      slotNzTotal[i]++;
      if (c.ex === 1) slotNzExact[i]++;
    }} else if (c.ex === 1) {{
      slotZeroExact[i]++;
    }}
  }});

  let thead = `<thead><tr><th class="row-hdr">Group</th>`;
  for (let s = 1; s <= 5; s++) {{
    thead += `<th>Slot ${{s}}</th>`;
    if (showSpacers && s === spacerAfterSlot) {{
      thead += `<th class="spacer-hdr">Unaligned<br/>0-Spacer</th>`;
    }}
  }}
  thead += `</tr></thead>`;

  let tbody = "<tbody>";
  for (let g = 1; g <= 21; g++) {{
    const khipuAG = activeTab === "part1" ? g : g + 24;
    const label = activeTab === "part1" ? `G${{g}}` : `A:G${{khipuAG}} / #${{g}}`;
    tbody += `<tr><th class="row-hdr">${{label}}</th>`;
    for (let s = 1; s <= 5; s++) {{
      const key = `${{g}}-${{s}}`;
      const c = byKey.get(key);
      if (!c) {{
        // Explicit unaligned G1 S1 (absent in Khipu A; unaligned p1 in Khipu B)
        tbody += `<td class="cell unaligned" title="Absent in Khipu A; Khipu B cord p1 is unaligned (excluded from denominator)"><span class="c-main">Unaligned</span><span class="c-sub">A:— · B:p1</span></td>`;
      }} else {{
        const cls = cellClass(c);
        const dim = passesFilter(c) ? "" : " dim";
        const sel = selectedCellKey === key ? " selected" : "";
        const pat = (activeTab === "part1" && c.s === 4 && c.g >= 1 && c.g <= 6) ? " pattern-ring" : "";
        const aria = `Group ${{g}} Slot ${{s}}: Khipu A value ${{c.av}} vs ${{activeTab === "part1" ? "Khipu B" : "Khipu C"}} value ${{c.bv}}, diff ${{c.d}}`;
        tbody += `<td class="cell ${{cls}}${{dim}}${{sel}}${{pat}}" data-key="${{key}}" tabindex="0" role="button" aria-label="${{aria}}">${{cellHtml(c)}}</td>`;
      }}
      if (showSpacers && s === spacerAfterSlot) {{
        const sp = spacers[String(g)];
        if (sp !== undefined) {{
          tbody += `<td class="cell unaligned" title="Unaligned structural 0-spacer cord p${{sp}} in ${{activeTab === 'part1' ? 'Khipu B' : 'Khipu C'}}"><span class="c-main">p${{sp}} = 0</span><span class="c-sub">spacer</span></td>`;
        }} else {{
          tbody += `<td class="cell unaligned" style="opacity:0.45;"><span class="c-main">—</span><span class="c-sub">no spacer</span></td>`;
        }}
      }}
    }}
    tbody += "</tr>";
  }}
  tbody += "</tbody>";

  let tfoot = `<tfoot>`;
  tfoot += `<tr class="nz-row"><td style="text-align:right;padding-right:7px;">Non-zero exact</td>`;
  for (let s = 1; s <= 5; s++) {{
    const i = s - 1;
    const pct = slotNzTotal[i] ? ((100 * slotNzExact[i]) / slotNzTotal[i]).toFixed(0) : "—";
    tfoot += `<td>${{slotNzExact[i]}}/${{slotNzTotal[i]}} (${{pct}}%)</td>`;
    if (showSpacers && s === spacerAfterSlot) tfoot += `<td>unaligned</td>`;
  }}
  tfoot += `</tr><tr><td style="text-align:right;padding-right:7px;">Zero exact (0=0)</td>`;
  for (let s = 1; s <= 5; s++) {{
    const i = s - 1;
    tfoot += `<td>${{slotZeroExact[i]}} zero(s)</td>`;
    if (showSpacers && s === spacerAfterSlot) tfoot += `<td>${{Object.keys(spacers).length}} spacers</td>`;
  }}
  tfoot += `</tr><tr><td style="text-align:right;padding-right:7px;">All aligned</td>`;
  for (let s = 1; s <= 5; s++) {{
    const i = s - 1;
    const pct = ((100 * slotAllExact[i]) / slotTotal[i]).toFixed(0);
    tfoot += `<td>${{slotAllExact[i]}}/${{slotTotal[i]}} (${{pct}}%)</td>`;
    if (showSpacers && s === spacerAfterSlot) tfoot += `<td>excl.</td>`;
  }}
  tfoot += `</tr></tfoot>`;

  const tbl = document.getElementById("matrixTable");
  tbl.innerHTML = thead + tbody + tfoot;

  tbl.querySelectorAll("td.cell[data-key]").forEach(td => {{
    const key = td.dataset.key;
    const c = byKey.get(key);
    const activate = () => {{
      selectedCellKey = selectedCellKey === key ? null : key;
      lastInspectedKey = key;
      renderMatrix();
      renderInspector(c);
    }};
    td.addEventListener("mouseenter", () => {{
      if (!selectedCellKey) {{
        lastInspectedKey = key;
        renderInspector(c);
      }}
    }});
    td.addEventListener("click", activate);
    td.addEventListener("keydown", (e) => {{
      if (e.key === "Enter" || e.key === " ") {{
        e.preventDefault();
        activate();
      }}
    }});
  }});

  // Preserve previously inspected cell across filter/mode changes instead of jumping to rows[0]
  const keepKey = selectedCellKey || lastInspectedKey;
  if (keepKey && byKey.has(keepKey)) {{
    renderInspector(byKey.get(keepKey));
  }} else if (rows.length > 0) {{
    lastInspectedKey = `${{rows[0].g}}-${{rows[0].s}}`;
    renderInspector(rows[0]);
  }}
}}

function renderPart2() {{
  const tbl = document.getElementById("part2Table");
  let html = `<thead><tr>
    <th>Check</th>
    <th>Target</th>
    <th>Part 2 Values</th>
    <th>Comparison / Best Window</th>
    <th class="num">Exact</th>
    <th class="num">Pearson r</th>
    <th class="num">Spearman ρ</th>
    <th class="num">Sum / Mean |Rel Err|</th>
    <th>Note</th>
  </tr></thead><tbody>`;
  DATA.part2.forEach(r => {{
    const badge = r.is_internal
      ? `<span class="pill" style="background:var(--exact-bg);color:var(--exact-ink);border-color:var(--exact-bd);">ρ = ${{r.spearman_rho.toFixed(3)}}</span>`
      : `<span class="pill warn">Absent (${{r.exact_matches}}/${{r.n_cells}})</span>`;
    const errTxt = r.is_internal
      ? `${{r.sum_rel_err_pct > 0 ? "+" : ""}}${{r.sum_rel_err_pct.toFixed(1)}}% / ${{r.mean_abs_rel_err_pct.toFixed(1)}}%`
      : "—";
    html += `<tr>
      <td><code>${{r.check}}</code></td>
      <td>${{r.target_khipu}}</td>
      <td><code>${{r.part2_vals}}</code></td>
      <td><code>${{r.comparison_vals}}</code></td>
      <td class="num">${{r.exact_matches}} / ${{r.n_cells}}</td>
      <td class="num">${{Number(r.pearson_r).toFixed(3)}}</td>
      <td class="num">${{badge}}</td>
      <td class="num">${{errTxt}}</td>
      <td>${{r.note}}</td>
    </tr>`;
  }});
  html += "</tbody>";
  tbl.innerHTML = html;
}}

function renderRevDup() {{
  const tbl = document.getElementById("revdupTable");
  let html = `<thead><tr>
    <th class="num">#</th>
    <th class="num">KH0121 Ord</th>
    <th class="num">KH0121 Val</th>
    <th class="num">KH0473 Ord (rev)</th>
    <th class="num">KH0473 Val (rev)</th>
    <th>Status</th>
    <th>Physical Length Parity (cm)</th>
    <th class="num">|Δlen|</th>
  </tr></thead><tbody>`;
  DATA.revdup.forEach(r => {{
    if (r.unaligned) {{
      html += `<tr class="unaligned-row">
        <td class="num">—</td>
        <td class="num">—</td>
        <td class="num"><em>absent</em></td>
        <td class="num"><strong>p${{r.o2}}</strong></td>
        <td class="num"><strong>0</strong></td>
        <td><span class="pill alert">Unaligned stub (museum check)</span></td>
        <td><span style="font-size:12px;">0.0 cm broken stub in KH0473 only (30 vs. 29 cords)</span></td>
        <td class="num">—</td>
      </tr>`;
      return;
    }}
    const isMajorDiff = r.ex === 0;
    const matchPill = r.ex === 1
      ? `<span class="pill" style="background:var(--exact-bg);color:var(--exact-ink);border-color:var(--exact-bd);">Exact</span>`
      : `<span class="pill alert">Discrepancy (museum check)</span>`;
    const w1 = Math.min(100, (r.l1 / 100) * 100).toFixed(1);
    const w2 = Math.min(100, (r.l2 / 100) * 100).toFixed(1);
    const barHtml = (r.l1 > 0 && r.l2 > 0) ? `
      <div class="mini-bars">
        <div class="mini-bar-row"><span style="width:42px;">${{r.l1.toFixed(1)}}</span><div class="mini-bar" style="width:${{w1}}%;background:#1d5b54;"></div></div>
        <div class="mini-bar-row"><span style="width:42px;">${{r.l2.toFixed(1)}}</span><div class="mini-bar" style="width:${{w2}}%;background:#c75d38;"></div></div>
      </div>` : `<span style="color:var(--muted);font-size:12px;">0.0 / 0.0 (unmeasured end)</span>`;
    html += `<tr class="${{isMajorDiff ? "museum-flag" : ""}}">
      <td class="num">${{r.idx}}</td>
      <td class="num">p${{r.o1}}</td>
      <td class="num"><strong>${{r.v1}}</strong></td>
      <td class="num">p${{r.o2}}</td>
      <td class="num"><strong>${{r.v2}}</strong></td>
      <td>${{matchPill}}</td>
      <td>${{barHtml}}</td>
      <td class="num">${{r.dl.toFixed(1)}} cm</td>
    </tr>`;
  }});
  html += "</tbody>";
  tbl.innerHTML = html;
}}

function switchTab(tab) {{
  activeTab = tab;
  selectedCellKey = null;
  lastInspectedKey = null;
  highlightSlot4 = false;
  document.getElementById("toggleSlot4Pattern").classList.remove("active");
  document.getElementById("toggleSlot4Pattern").setAttribute("aria-pressed", "false");
  document.querySelectorAll(".tab-btn").forEach(b => b.classList.toggle("active", b.dataset.tab === tab));
  document.getElementById("matrixView").classList.toggle("hidden", tab !== "part1" && tab !== "part3");
  document.getElementById("part2View").classList.toggle("hidden", tab !== "part2");
  document.getElementById("revdupView").classList.toggle("hidden", tab !== "revdup");
  if (tab === "part1" || tab === "part3") renderMatrix();
}}

document.querySelectorAll(".tab-btn").forEach(btn => {{
  btn.addEventListener("click", () => switchTab(btn.dataset.tab));
}});
document.querySelectorAll("#filterBtns .ctrl").forEach(btn => {{
  btn.addEventListener("click", () => {{
    highlightSlot4 = false;
    document.getElementById("toggleSlot4Pattern").classList.remove("active");
    document.getElementById("toggleSlot4Pattern").setAttribute("aria-pressed", "false");
    activeFilter = btn.dataset.filter;
    document.querySelectorAll("#filterBtns .ctrl").forEach(b => b.classList.toggle("active", b === btn));
    renderMatrix();
  }});
}});
document.querySelectorAll("#modeBtns .ctrl").forEach(btn => {{
  btn.addEventListener("click", () => {{
    activeMode = btn.dataset.mode;
    document.querySelectorAll("#modeBtns .ctrl").forEach(b => b.classList.toggle("active", b === btn));
    renderMatrix();
  }});
}});
document.getElementById("toggleSpacers").addEventListener("click", (e) => {{
  showSpacers = !showSpacers;
  e.currentTarget.classList.toggle("active", showSpacers);
  e.currentTarget.setAttribute("aria-pressed", String(showSpacers));
  renderMatrix();
}});
document.getElementById("toggleSlot4Pattern").addEventListener("click", (e) => {{
  highlightSlot4 = !highlightSlot4;
  e.currentTarget.classList.toggle("active", highlightSlot4);
  e.currentTarget.setAttribute("aria-pressed", String(highlightSlot4));
  renderMatrix();
}});

renderMatrix();
renderPart2();
renderRevDup();
</script>
</body>
</html>
"""


def main() -> None:
    data = build_data()
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(render_html(data), encoding="utf-8")
    print(f"Generated {OUT} ({OUT.stat().st_size:,} bytes)")
    print("Source CSV SHA-256 digests:")
    for k, v in data["checksums"].items():
        print(f"  {k}: {v}")
    print("Summaries:", data["summary"])


if __name__ == "__main__":
    main()
