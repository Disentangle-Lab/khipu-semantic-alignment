"""
make_figures.py — Publication-quality figures (300 DPI) for khipu-semantic-alignment:
  - fig1_step1_archive_invariance.png:   Step 1 (site, k) excess + matched n=2 subsampling CI
  - fig2_step2_taxonomy_outliers.png:    Step 2 (G, k) landscape, F1-F4 vs OUT, HDBSCAN noise (-1)
  - fig3_subsidiary_grammar_outliers.png: Exact intra-row subsidiary column grammar (76.5% BH-FDR)
                                          & deterministic outlier Fisher enrichment (OR = 7.4-10.7)
  - fig4_gate0_vs_step3_alignment.png:   Gate 0 positive-control calibration (A_JSD 0.977-0.985)
                                          vs pre-registered Step 3 null (0/30 & 0/22 FDR) + Santa Valley replication.
"""

from __future__ import annotations

import ast
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from src.data_loader import FIGURES_DIR, RESULTS_DIR

sns.set_theme(style="whitegrid", context="paper", font_scale=1.08)


def make_fig1():
    sub = pd.read_csv(RESULTS_DIR / "step1_matched_subsample.csv")
    tests = pd.read_csv(RESULTS_DIR / "step1_archive_tests.csv")
    if sub.empty:
        return

    fig, axes = plt.subplots(1, 2, figsize=(13.2, 5.2))
    for ax, metric, title, color in [
        (axes[0], "S1a_color_modal_orient_inv", "S1a: Column-Color Invariance (within same k)", "#2E86AB"),
        (axes[1], "S1b_num_rank_okr_orient_inv", "S1b: Column-Numeric Rank Invariance (within same k)", "#A23B72"),
    ]:
        d = sub[sub["metric"] == metric].copy()
        d["label"] = d.apply(lambda r: f"{r['site']} (k={r['k']}, n={r['n_cell']})", axis=1)
        d = d.sort_values("full_excess", ascending=True).reset_index(drop=True)
        y = np.arange(len(d))

        ax.axvline(0, color="#444444", linestyle="--", linewidth=1.0, alpha=0.8)
        lo_exc = d["sub2_ci_lo"] - d["between_k_mean"]
        hi_exc = d["sub2_ci_hi"] - d["between_k_mean"]
        med_exc = d["sub2_median"] - d["between_k_mean"]

        ax.hlines(y, lo_exc, hi_exc, color=color, alpha=0.45, linewidth=3.5, label="Matched n=2 (95% CI)")
        ax.scatter(med_exc, y, color=color, s=42, zorder=3, marker="o", label="Matched n=2 median")
        ax.scatter(d["full_excess"], y, color="#111111", s=34, zorder=4, marker="D", label="Full-cell mean")

        ax.set_yticks(y)
        ax.set_yticklabels(d["label"], fontsize=8.8)
        ax.set_xlabel("Within-archive minus between-archive similarity (same k)")

        trow = tests[tests["metric"] == f"primary_86__{metric}"]
        if not trow.empty:
            tm = trow.iloc[0]["T_macro_obs"]
            pm = trow.iloc[0]["p_macro"]
            ax.set_title(f"{title}\nT_macro = {tm:+.3f}, p = {pm:.4f}", fontweight="bold")
        else:
            ax.set_title(title, fontweight="bold")
        ax.legend(loc="lower right", fontsize=8, frameon=True)

    plt.tight_layout()
    out = FIGURES_DIR / "fig1_step1_archive_invariance.png"
    plt.savefig(out, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Saved {out}")


def make_fig2():
    sig = pd.read_csv(RESULTS_DIR / "step2_signatures_86.csv")
    if sig.empty:
        return

    fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.5))

    fam_palette = {
        "F1_compact_parallel": "#2E86AB",
        "F2_medium_multi_attr": "#28A745",
        "F3_self_totaling": "#F18F01",
        "F4_nonclean_hybrid": "#8E44AD",
        "OUT": "#D62828",
    }
    ax = axes[0]
    for fam in ["F1_compact_parallel", "F2_medium_multi_attr", "F3_self_totaling", "F4_nonclean_hybrid", "OUT"]:
        grp = sig[sig["det_family"] == fam]
        if grp.empty:
            continue
        ax.scatter(
            grp["k"], grp["G"],
            s=36 + 115 * grp["sub_ratio"].clip(0, 1.5),
            alpha=0.84,
            color=fam_palette[fam],
            edgecolors="black",
            linewidths=0.6 if fam != "OUT" else 1.1,
            label=f"{fam} (n={len(grp)})",
            zorder=4 if fam == "OUT" else 3,
        )

    # Clean, non-overlapping annotations for key landmark outliers
    landmark_offsets = {
        "KH0128": ("KH0128 (2×160)", (-72, 18)),
        "KH0286": ("KH0286 (2×36)", (-10, 22)),
        "KH0242": ("KH0242 Calendar (13×22)", (12, 10)),
        "KH0300": ("KH0300 (3×20)", (14, 6)),
        "KH0097": ("KH0097 (10×14)", (16, -14)),
        "KH0237": ("KH0237 (6×16)", (16, 4)),
        "KH0183": ("KH0183 (19×11)", (12, 10)),
    }
    for okr, (lbl, offset) in landmark_offsets.items():
        row = sig[sig["OKR_NUM"] == okr]
        if not row.empty:
            x, y = float(row.iloc[0]["k"]), float(row.iloc[0]["G"])
            ax.annotate(
                lbl,
                xy=(x, y),
                xytext=offset,
                textcoords="offset points",
                fontsize=7.8,
                color="#780000",
                fontweight="bold" if okr in ("KH0242", "KH0128") else "normal",
                arrowprops=dict(arrowstyle="-", color="#999999", lw=0.7),
            )

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(2.5, 220)
    ax.set_ylim(1.7, 90)
    ax.set_xlabel("Columns per modal group (k = h3_m, log scale)")
    ax.set_ylabel("Number of modal groups (G = h3_G, log scale)")
    ax.set_title("86 Doubly-Tabular Khipus: Deterministic Families & Outliers\n(point size = subsidiary ratio)", fontweight="bold")
    ax.legend(fontsize=8, loc="upper right", frameon=True)

    # Right panel: H3a vs H3b with outlier overlap
    ax2 = axes[1]
    inliers = sig[(~sig["is_det_outlier"]) & (~sig["is_hdbscan_noise"])]
    det_only = sig[(sig["is_det_outlier"]) & (~sig["is_hdbscan_noise"])]
    hdb_only = sig[(~sig["is_det_outlier"]) & (sig["is_hdbscan_noise"])]
    both_out = sig[(sig["is_det_outlier"]) & (sig["is_hdbscan_noise"])]

    ax2.scatter(inliers["h3_col_match"], inliers["h3_val_r2"], s=45, alpha=0.68, color="#3A7CA5", label=f"Inlier both (n={len(inliers)})")
    ax2.scatter(hdb_only["h3_col_match"], hdb_only["h3_val_r2"], s=55, alpha=0.85, color="#F18F01", marker="^", edgecolors="black", linewidths=0.5, label=f"HDBSCAN noise only (n={len(hdb_only)})")
    ax2.scatter(det_only["h3_col_match"], det_only["h3_val_r2"], s=60, alpha=0.85, color="#D62828", marker="s", edgecolors="black", linewidths=0.5, label=f"Det OUT only (n={len(det_only)})")
    ax2.scatter(both_out["h3_col_match"], both_out["h3_val_r2"], s=78, alpha=0.95, color="#780000", marker="X", edgecolors="black", linewidths=0.6, label=f"Both OUT & HDBSCAN -1 (n={len(both_out)})")

    ax2.set_xlabel("H3a Color Positional Homology (match rate)")
    ax2.set_ylabel("H3b Numeric Positional Homology (R²)")
    ax2.set_title("Outlier Accounting Across Two Independent Criteria\n(Det OUT: 24.4% | HDBSCAN noise: 52.3%)", fontweight="bold")
    ax2.legend(fontsize=8, loc="lower right", frameon=True)

    plt.tight_layout()
    out = FIGURES_DIR / "fig2_step2_taxonomy_outliers.png"
    plt.savefig(out, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Saved {out}")


def make_fig3():
    sub_perm = pd.read_csv(RESULTS_DIR / "step2_subsidiary_perm_5000.csv")
    enrich = pd.read_csv(RESULTS_DIR / "step2_outlier_enrichment_fisher.csv")

    fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.2))

    # Left panel: Exact intra-row permutation subsidiary column attachment (binary presence, n=34)
    ax = axes[0]
    b = sub_perm[sub_perm["cords_with_sub"] >= 5].copy().sort_values("ent_binary", ascending=True).reset_index(drop=True)
    colors = ["#1B998B" if q < 0.05 else "#BBBBBB" for q in b["q_perm_binary_var"]]
    y = np.arange(len(b))
    ax.barh(y, b["ent_binary"], color=colors, edgecolor="none", height=0.78)
    for yi, row in zip(y, b.itertuples(index=False)):
        if row.ent_binary == 0.0:
            ax.scatter([0.014], [yi], color="#1B998B", s=26, zorder=4)
            ax.text(0.032, yi, "0.00 (single slot)", va="center", ha="left", fontsize=7.4, color="#136F63", fontweight="bold")
    ax.axvline(b["ent_binary"].median(), color="#111111", linestyle="--", lw=1.1, label=f"Median entropy = {b['ent_binary'].median():.3f}")
    ax.axvline(1.0, color="#888888", linestyle=":", lw=1.0, label="Uniform across columns (1.0)")

    # Highlight iconic specimens on y-axis ticks (keep every label clean)
    labels = [
        f"{r.OKR_NUM} ({r.G}×{r.k})" if r.OKR_NUM in {"KH0156", "KH0469", "KH0192", "KH0068", "KH0242", "KH0111", "KH0287", "KH0247"} else ""
        for r in b.itertuples(index=False)
    ]
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=8.0)
    ax.set_xlim(0, 1.08)
    ax.set_xlabel("Normalized Column Slot Entropy of Subsidiary Presence (0 = single slot, 1 = uniform)")
    n_sig = int((b["q_perm_binary_var"] < 0.05).sum())
    ax.set_title(
        f"Subsidiary Column Grammar (Exact Intra-Row Perm., N=5,000)\n"
        f"BH-FDR q < 0.05: {n_sig} / {len(b)} ({100*n_sig/len(b):.1f}%, binom p = $2.84 \\times 10^{{-29}}$)",
        fontweight="bold",
    )
    from matplotlib.patches import Patch
    handles, leg_labels = ax.get_legend_handles_labels()
    handles = [Patch(facecolor="#1B998B", label=f"BH q < 0.05 (n={n_sig})"), Patch(facecolor="#BBBBBB", label=f"Not sig. (n={len(b)-n_sig})")] + handles
    ax.legend(handles=handles, fontsize=8, loc="lower right", frameon=True)

    # Right panel: Fisher exact enrichment of Deterministic OUT vs HDBSCAN noise
    ax2 = axes[1]
    targets = [
        ("Official calendar + UR006/KH0242", "Official Calendars\n(KFG10 + UR006)"),
        ("KFG Part B new paired khipus (9 pairs)", "New Paired Khipus\n(KFG Part B)"),
        ("UNION (any calendar OR new pair)", "Union\n(Calendar or Pair)"),
    ]
    x = np.arange(len(targets))
    width = 0.34
    det_ors, det_ps, hdb_ors, hdb_ps = [], [], [], []
    for t_key, _ in targets:
        rd = enrich[(enrich["outlier_def"] == "Deterministic OUT (21/86)") & (enrich["target"] == t_key)].iloc[0]
        rh = enrich[(enrich["outlier_def"] == "HDBSCAN noise -1 (45/86)") & (enrich["target"] == t_key)].iloc[0]
        det_ors.append(rd["odds_ratio"])
        det_ps.append(rd["fisher_p_two"])
        hdb_ors.append(rh["odds_ratio"])
        hdb_ps.append(rh["fisher_p_two"])

    b1 = ax2.bar(x - width / 2, det_ors, width, color="#D62828", alpha=0.88, edgecolor="black", lw=0.6, label="Deterministic OUT (21/86)")
    b2 = ax2.bar(x + width / 2, hdb_ors, width, color="#F18F01", alpha=0.78, edgecolor="black", lw=0.6, label="HDBSCAN noise -1 (45/86)")
    ax2.axhline(1.0, color="#333333", linestyle="--", lw=1.0, label="No enrichment (OR = 1.0)")

    for rect, orv, pv in zip(b1, det_ors, det_ps):
        star = "*" if pv < 0.05 else ""
        ax2.text(rect.get_x() + rect.get_width() / 2, orv + 0.35, f"OR={orv:.1f}{star}\np={pv:.3f}", ha="center", va="bottom", fontsize=8.2, fontweight="bold" if pv < 0.05 else "normal", color="#780000")
    for rect, orv, pv in zip(b2, hdb_ors, hdb_ps):
        y_txt = max(orv + 0.35, 1.22)
        ax2.text(
            rect.get_x() + rect.get_width() / 2, y_txt,
            f"OR={orv:.2f}\np={pv:.2f}",
            ha="center", va="bottom", fontsize=8.0, color="#7A4900",
            bbox=dict(facecolor="white", edgecolor="none", pad=0.6, alpha=0.85),
        )

    ax2.set_xticks(x)
    ax2.set_xticklabels([t[1] for t in targets], fontsize=9.2)
    ax2.set_ylim(0, max(det_ors) * 1.28)
    ax2.set_ylabel("Odds Ratio (Outliers vs. Inliers)")
    ax2.set_title("Functional Enrichment of Outliers (Fisher Exact)\nGeometric Outliers Capture Calendars & Pairs; HDBSCAN Noise Does Not", fontweight="bold")
    ax2.legend(fontsize=8.2, loc="upper right", frameon=True)

    plt.tight_layout()
    out = FIGURES_DIR / "fig3_subsidiary_grammar_outliers.png"
    plt.savefig(out, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Saved {out}")


def make_fig4():
    g0a = pd.read_csv(RESULTS_DIR / "gate0a_santa_valley_replication.csv")
    g0b = pd.read_csv(RESULTS_DIR / "gate0b_santa_valley_calibration.csv")
    s3 = pd.read_csv(RESULTS_DIR / "step3_alignment_pairs.csv")

    fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.2))

    # Left panel: Positive Controls (Gate 0B + Text-to-Text Huaraz) vs Pre-Registered Step 3 (all schemas)
    ax = axes[0]
    pos_vals = g0b[g0b["khipu_b"] != "UR092"]["A_JSD"].tolist() + [0.9774]  # 4 Santa Valley core pairs + Huaraz 1532 vs 1558
    prim_all = s3[s3["cohort"] == "primary_86"]["A_JSD"].tolist()
    sens_all = s3[s3["cohort"] == "sens_H3b_only"]["A_JSD"].tolist()

    groups = [
        ("Positive Controls\n(Santa Valley k=6 &\nHuaraz 1532 vs. 1558)", pos_vals, "#1B998B"),
        (f"Step 3 Primary 86\nvs. Textos Andinos\n(k=5,6,7,22; 0/{len(prim_all)} FDR)", prim_all, "#3A7CA5"),
        (f"Step 3 H3b-only\nvs. Textos Andinos\n(k=5,6,7; 0/{len(sens_all)} FDR)", sens_all, "#8E44AD"),
    ]
    rng = np.random.default_rng(42)
    for i, (lbl, vals, col) in enumerate(groups):
        vals = np.asarray(vals)
        bp = ax.boxplot(
            vals, positions=[i], widths=0.42, patch_artist=True, showfliers=False,
            boxprops=dict(facecolor=col, alpha=0.28, edgecolor=col, lw=1.2),
            medianprops=dict(color="#111111", lw=1.8),
            whiskerprops=dict(color=col, lw=1.1),
            capprops=dict(color=col, lw=1.1),
        )
        jit = rng.uniform(-0.11, 0.11, size=len(vals))
        ax.scatter(np.full(len(vals), i) + jit, vals, color=col, s=36, alpha=0.88, edgecolors="black", linewidths=0.45, zorder=4)

    ax.axhline(0.977, color="#1B998B", linestyle="--", lw=1.0, alpha=0.75, label=r"Positive-control floor ($A_{\mathrm{JSD}} = 0.977$)")
    ax.set_xticks(range(len(groups)))
    ax.set_xticklabels([g[0] for g in groups], fontsize=8.8)
    ymin = min(0.50, min(prim_all + sens_all) - 0.04)
    ax.set_ylim(ymin, 1.015)
    ax.set_ylabel(r"Pre-registered Alignment Score $A_{\mathrm{JSD}} = 1 - \mathrm{JSD}_2$")
    ax.set_title("Gate 0 Positive Controls vs. Pre-Registered Step 3\nMetric Detects Same-Archive/Region Pairs; Cross-Region Null (0 FDR)", fontweight="bold")
    ax.legend(fontsize=8.2, loc="lower left", frameon=True)

    # Right panel: Gate 0A Santa Valley 1670 Exact Replication (First-cord tribute means by moiety)
    ax2 = axes[1]
    order = ["UR092", "UR087", "UR091", "UR089", "UR090", "UR088"]
    g0a_ord = g0a.set_index("kfg_name").loc[order].reset_index()
    moiety_colors = ["#2E86AB", "#2E86AB", "#D62828", "#D62828", "#F18F01", "#D62828"]
    x2 = np.arange(len(g0a_ord))
    bars = ax2.bar(x2, g0a_ord["first_cord_mean_kfg"], color=moiety_colors, edgecolor="black", lw=0.6, width=0.62, alpha=0.88)

    for rect, row in zip(bars, g0a_ord.itertuples(index=False)):
        ax2.text(
            rect.get_x() + rect.get_width() / 2,
            rect.get_height() + 0.22,
            f"{row.first_cord_mean_kfg:.2f}\n(n={row.n_groups_obs})",
            ha="center", va="bottom", fontsize=8.0,
        )

    from matplotlib.patches import Patch
    ax2.legend(
        handles=[
            Patch(facecolor="#2E86AB", edgecolor="black", lw=0.6, label="Recto moiety (UR092, UR087)"),
            Patch(facecolor="#D62828", edgecolor="black", lw=0.6, label="Verso moiety (UR091, UR089, UR088)"),
            Patch(facecolor="#F18F01", edgecolor="black", lw=0.6, label="Mixed R/V (UR090)"),
        ],
        fontsize=8.0, loc="upper left", frameon=True,
    )
    ax2.set_xticks(x2)
    ax2.set_xticklabels([f"{r.kfg_name}\n({r.OKR_NUM})" for r in g0a_ord.itertuples(index=False)], fontsize=8.5)
    ax2.set_ylim(0, 10.6)
    ax2.set_ylabel("Mean First-Cord Tribute Value (pesos, KFG Excel)")
    ax2.set_title(
        "Gate 0A: Santa Valley 1670 Computational Replication\n"
        "133/133 Groups | 63R/69V | First-Cord Sum = 369 vs. 367 Pesos",
        fontweight="bold",
    )

    plt.tight_layout()
    out = FIGURES_DIR / "fig4_gate0_vs_step3_alignment.png"
    plt.savefig(out, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Saved {out}")


def main():
    make_fig1()
    make_fig2()
    make_fig3()
    make_fig4()


if __name__ == "__main__":
    main()
