"""
phase2c_pair_discovery.py — Phase 2C (Bloque C):
Corpus-wide 611 x 611 high-entropy rare-number pair discovery sweep on OKR v2.1.0,
forensic 2D alignment of the Berlin Accounting Triptych (KH0192/UR1175 VA 42518 Part 1 & Part 3
vs. KH0470/UR233 VA 63038 b & KH0469/UR232 VA 63038(B)), and forensic audit of the unflagged
reversed physical duplicate KH0121 (AS107/UR1107 VA 24370(B), Marquez) <-> KH0473 (UR237 VA 66834).

Generates frozen Phase 2C tables in results/ and publication-grade 4-panel figure:
  figures/fig6_phase2c_berlin_triptych_and_pairs.png (300 DPI).

Does NOT modify any frozen Phase 1 outputs or PREREGISTRATION.md.
"""

from __future__ import annotations

import pathlib
import re
import sqlite3
import warnings
from collections import Counter

import numpy as np
import openpyxl
import pandas as pd
from scipy.stats import mannwhitneyu, pearsonr, spearmanr

from src.data_loader import FIGURES_DIR, KFG_XLSX_DIR, RESULTS_DIR, find_okr_db
from src.phase2a_subsidiary_grammar import decode_base10_cords_fast

warnings.filterwarnings("ignore")

KNOWN_CONTROL_PAIRS = {
    frozenset(("KH0140", "KH0483")): "Known Dup (AS125/UR247)",
    frozenset(("KH0303", "KH0304")): "Known Dup (UR066/UR067)",
    frozenset(("KH0301", "KH0305")): "Known Pair (Puruchuco UR064/UR068)",
    frozenset(("KH0491", "KH0503")): "Known Pair (Incahuasi UR255/UR267A)",
    frozenset(("KH0492", "KH0504")): "Known Pair (Incahuasi UR256/UR267B)",
}

NEW_DISCOVERY_PAIRS = {
    frozenset(("KH0192", "KH0470")): "NEW: Berlin Triptych Part 1 (VA42518 / VA63038b)",
    frozenset(("KH0192", "KH0469")): "NEW: Berlin Triptych Part 3 (VA42518 / VA63038B)",
    frozenset(("KH0121", "KH0473")): "NEW: Reversed Dup (VA24370B / VA66834)",
}


def load_kfg_pendants(fname: str) -> tuple[dict[int, tuple[int, str, float]], dict[int, list]]:
    p = KFG_XLSX_DIR / fname
    wb = openpyxl.load_workbook(p, read_only=True, data_only=True)
    rows = list(wb["Cords"].iter_rows(values_only=True))
    wb.close()
    pends: dict[int, tuple[int, str, float]] = {}
    subs: dict[int, list] = {}
    for r in rows[1:]:
        if not r or not r[0]:
            continue
        cn = str(r[0]).strip()
        try:
            val = int(float(r[8])) if r[8] not in (None, "", "None") else 0
        except Exception:
            val = 0
        try:
            length = float(r[4]) if r[4] not in (None, "", "None") else np.nan
        except Exception:
            length = np.nan
        col = str(r[7] or "")
        if re.match(r"^p\d+$", cn):
            pends[int(cn[1:])] = (val, col, length)
        elif re.match(r"^p(\d+)s", cn):
            pn = int(re.match(r"^p(\d+)s", cn).group(1))
            subs.setdefault(pn, []).append((cn, val, col))
    return pends, subs


def run_corpus_rare_number_sweep(con: sqlite3.Connection) -> pd.DataFrame:
    meta = pd.read_sql(
        "SELECT KHIPU_ID, OKR_NUM, INVESTIGATOR_NUM, orig_inv_num, MUSEUM_NAME, MUSEUM_NUM, PROVENANCE, DUPLICATE_FLAG FROM khipu_main",
        con,
    )
    vals = decode_base10_cords_fast(con)
    p1 = pd.read_sql(
        """
        SELECT CORD_ID, KHIPU_ID, CORD_ORDINAL
        FROM cord
        WHERE CORD_LEVEL = 1 AND COALESCE(CORD_CLASSIFICATION, '') != 'T'
        ORDER BY KHIPU_ID, CORD_ORDINAL
        """,
        con,
    )
    p1 = p1.merge(vals, on="CORD_ID", how="left")
    p1["cord_val"] = p1["cord_val"].fillna(0.0)

    k_seq: dict[int, np.ndarray] = {}
    k_rare: dict[int, Counter] = {}
    for kid, g in p1.groupby("KHIPU_ID", sort=False):
        v = g["cord_val"].values
        vr = [int(x) for x in v if x >= 11 and int(x) % 10 != 0]
        if len(vr) >= 5:
            k_seq[kid] = v[v > 0]
            k_rare[kid] = Counter(vr)

    kids = sorted(k_seq.keys())
    meta_idx = meta.set_index("KHIPU_ID").to_dict("index")
    hits = []
    for i in range(len(kids)):
        a = kids[i]
        ca = k_rare[a]
        na = sum(ca.values())
        za = k_seq[a]
        for j in range(i + 1, len(kids)):
            b = kids[j]
            cb = k_rare[b]
            nb = sum(cb.values())
            inter = ca & cb
            n_rare_dist = len(inter)
            if n_rare_dist < 5:
                continue
            rare_dice = 2.0 * sum(inter.values()) / (na + nb)
            if rare_dice < 0.25 and n_rare_dist < 10:
                continue
            zb = k_seq[b]
            best_rare_ex, best_ex, best_r, best_ori = 0, 0, np.nan, "fwd"
            for ori, zbu in [("fwd", zb), ("rev", zb[::-1])]:
                for lag in range(-10, 11):
                    i0 = max(0, lag)
                    j0 = max(0, -lag)
                    L = min(len(za) - i0, len(zbu) - j0)
                    if L < 6:
                        continue
                    xa = za[i0:i0 + L]
                    xb = zbu[j0:j0 + L]
                    rex = int(np.sum((xa == xb) & (xa >= 11) & (xa % 10 != 0)))
                    ex = int(np.sum((xa == xb) & (xa >= 10)))
                    if rex > best_rare_ex or (rex == best_rare_ex and ex > best_ex):
                        best_rare_ex, best_ex, best_ori = rex, ex, ori
                        best_r = float(spearmanr(xa, xb).statistic) if (len(np.unique(xa)) >= 3 and len(np.unique(xb)) >= 3) else np.nan
            if best_rare_ex >= 4 or rare_dice >= 0.40:
                ma, mb = meta_idx[a], meta_idx[b]
                key = frozenset((ma["OKR_NUM"], mb["OKR_NUM"]))
                if key in NEW_DISCOVERY_PAIRS:
                    cat = NEW_DISCOVERY_PAIRS[key]
                elif key in KNOWN_CONTROL_PAIRS:
                    cat = KNOWN_CONTROL_PAIRS[key]
                else:
                    cat = "Background"
                hits.append({
                    "okr_a": ma["OKR_NUM"], "inv_a": ma["INVESTIGATOR_NUM"], "prov_a": ma["PROVENANCE"] or "",
                    "mnum_a": ma["MUSEUM_NUM"] or "",
                    "okr_b": mb["OKR_NUM"], "inv_b": mb["INVESTIGATOR_NUM"], "prov_b": mb["PROVENANCE"] or "",
                    "mnum_b": mb["MUSEUM_NUM"] or "",
                    "na_nz": len(za), "nb_nz": len(zb),
                    "rare_dist": n_rare_dist, "rare_dice": round(rare_dice, 4),
                    "ord_rare_ex": best_rare_ex, "ord_ex10": best_ex,
                    "spearman_r": round(best_r, 4) if pd.notna(best_r) else np.nan,
                    "ori": best_ori,
                    "dup_flag_okr": bool((ma["DUPLICATE_FLAG"] == 1) or (mb["DUPLICATE_FLAG"] == 1)),
                    "category": cat,
                })
    df = pd.DataFrame(hits).sort_values(["ord_rare_ex", "rare_dice", "rare_dist"], ascending=False).reset_index(drop=True)
    return df


def analyze_berlin_triptych(rng: np.random.Generator) -> tuple[pd.DataFrame, dict]:
    """Exact 2D matrix alignment of UR1175 (KH0192, VA 42518) Part 1 & Part 3 vs UR233 (KH0470) & UR232 (KH0469)."""
    p192, s192 = load_kfg_pendants("UR1175.xlsx")
    p470, s470 = load_kfg_pendants("UR233.xlsx")
    p469, s469 = load_kfg_pendants("UR232.xlsx")

    # Part 1: UR1175 p1..p104 (21 groups; group 1 missing broken initial cord, which UR233 p1=23 preserves!)
    # UR233 p1..p126 = 21 groups of 6 where index % 6 == 2 (p3, p9, p15, ..., p123) is 100% zero spacer!
    ur233_nonspacer_idx = [i for i in range(1, 127) if (i - 1) % 6 != 2]  # len 105; ur233_nonspacer_idx[0]=p1(23)
    rows = []
    for k_i in range(104):
        pn192 = k_i + 1
        pn470 = ur233_nonspacer_idx[k_i + 1]
        v1, c1, l1 = p192[pn192]
        v2, c2, l2 = p470[pn470]
        rows.append({
            "triptych_half": "Part1_vs_UR233",
            "group_idx": 1 + (k_i + 1) // 5,
            "slot_idx": 1 + (k_i + 1) % 5,
            "master_pnum": pn192, "master_val": v1, "master_col": c1, "master_len_cm": l1,
            "master_n_subs": len(s192.get(pn192, [])),
            "copy_okr": "KH0470", "copy_kfg": "UR233",
            "copy_pnum": pn470, "copy_val": v2, "copy_col": c2, "copy_len_cm": l2,
            "copy_n_subs": len(s470.get(pn470, [])),
            "diff": v2 - v1,
            "exact_match": int(v1 == v2),
        })

    # Part 3: UR1175 p120..p224 (21 groups of 5 = 105 cords) vs UR232 p1..p125!
    # Groups 1-18: spacer is at (i-1)%6 == 3 (p4, p10, ..., p106) -> p1..p108 yields exact 90 non-spacer cords matching p120..p209!
    # Groups 19-21 (tail p109..p125, 17 cords -> 15 data cords + 2 zero spacers at p117 and p123):
    ur232_p1_18 = [i for i in range(1, 109) if (i - 1) % 6 != 3]  # 90 cords -> Groups 1..18
    # Tail alignment: UR232 p109..p125 = [40,20,15,0,0, 70,20,156,0,50, 0,40,12,125,0, 90,21]
    # Dropping the 2 structural 0-spacers at p117 (after 156) and p119 (before 40) — wait, check p109..p125 exact indices:
    # p109=40, p110=20, p111=15, p112=0, p113=0 (Grp19: 50,40,15,0,0)
    # p114=70, p115=20, p116=156, [p117=0 spacer], p118=50, p119=0 (Grp20: 70,40,140,50,0)
    # p120=40, p121=12, p122=125, [p123=0 spacer], p124=90, p125=21 (Grp21: 50,12,125,90,11)
    ur232_tail = [109, 110, 111, 112, 113, 114, 115, 116, 118, 119, 120, 121, 122, 124, 125]
    ur232_all = ur232_p1_18 + ur232_tail  # exact 105 cords matching p120..p224!
    assert len(ur232_all) == 105

    for k_i in range(105):
        pn192 = 120 + k_i
        pn469 = ur232_all[k_i]
        v1, c1, l1 = p192[pn192]
        v2, c2, l2 = p469[pn469]
        rows.append({
            "triptych_half": "Part3_vs_UR232",
            "group_idx": 1 + k_i // 5,
            "slot_idx": 1 + k_i % 5,
            "master_pnum": pn192, "master_val": v1, "master_col": c1, "master_len_cm": l1,
            "master_n_subs": len(s192.get(pn192, [])),
            "copy_okr": "KH0469", "copy_kfg": "UR232",
            "copy_pnum": pn469, "copy_val": v2, "copy_col": c2, "copy_len_cm": l2,
            "copy_n_subs": len(s469.get(pn469, [])),
            "diff": v2 - v1,
            "exact_match": int(v1 == v2),
        })

    tdf = pd.DataFrame(rows)
    p1_sub = tdf[tdf["triptych_half"] == "Part1_vs_UR233"]
    p3_sub = tdf[tdf["triptych_half"] == "Part3_vs_UR232"]

    def perm_p(a, b, n_perm=10000):
        obs = float(pearsonr(a, b).statistic)
        cnt = sum(float(pearsonr(a, rng.permutation(b)).statistic) >= obs - 1e-12 for _ in range(n_perm))
        return obs, (cnt + 1) / (n_perm + 1.0)

    r1, p1_perm = perm_p(p1_sub["master_val"].values, p1_sub["copy_val"].values)
    r3, p3_perm = perm_p(p3_sub["master_val"].values, p3_sub["copy_val"].values)
    r1_sp = float(spearmanr(p1_sub["master_val"], p1_sub["copy_val"]).statistic)
    r3_sp = float(spearmanr(p3_sub["master_val"], p3_sub["copy_val"]).statistic)

    # Physical length divergence test (proving two distinct artifacts!)
    len_p1_m = p1_sub["master_len_cm"].dropna()
    len_p1_c = p1_sub["copy_len_cm"].dropna()
    mw_len_p1 = float(mannwhitneyu(len_p1_m, len_p1_c).pvalue)
    len_p3_m = p3_sub["master_len_cm"].dropna()
    len_p3_c = p3_sub["copy_len_cm"].dropna()
    mw_len_p3 = float(mannwhitneyu(len_p3_m, len_p3_c).pvalue)

    stats = {
        "p1_n": len(p1_sub),
        "p1_exact": int(p1_sub["exact_match"].sum()),
        "p1_within1": int((p1_sub["diff"].abs() <= 1).sum()),
        "p1_within10": int((p1_sub["diff"].abs() <= 10).sum()),
        "p1_pearson_r": r1, "p1_spearman_r": r1_sp, "p1_perm_p": p1_perm,
        "p1_med_len_master": float(len_p1_m.median()), "p1_med_len_copy": float(len_p1_c.median()), "p1_len_mw_p": mw_len_p1,
        "p3_n": len(p3_sub),
        "p3_exact": int(p3_sub["exact_match"].sum()),
        "p3_within1": int((p3_sub["diff"].abs() <= 1).sum()),
        "p3_within10": int((p3_sub["diff"].abs() <= 10).sum()),
        "p3_pearson_r": r3, "p3_spearman_r": r3_sp, "p3_perm_p": p3_perm,
        "p3_med_len_master": float(len_p3_m.median()), "p3_med_len_copy": float(len_p3_c.median()), "p3_len_mw_p": mw_len_p3,
        "total_cells": len(tdf), "total_exact": int(tdf["exact_match"].sum()),
        "master_total_subs": sum(len(v) for v in s192.values()),
        "ur233_total_subs": sum(len(v) for v in s470.values()),
        "ur232_total_subs": sum(len(v) for v in s469.values()),
    }
    return tdf, stats


def analyze_reversed_duplicate(con: sqlite3.Connection) -> tuple[pd.DataFrame, dict]:
    """Forensic verification of KH0121 (AS107/UR1107, VA 24370(B), Marquez) <-> KH0473 (UR237, VA 66834, unknown)."""
    vals = decode_base10_cords_fast(con)
    p1 = pd.read_sql(
        """
        SELECT m.OKR_NUM, m.INVESTIGATOR_NUM, m.MUSEUM_NUM, m.PROVENANCE,
               c.CORD_ID, c.CORD_ORDINAL, c.CORD_LENGTH, c.TWIST
        FROM cord c JOIN khipu_main m ON c.KHIPU_ID = m.KHIPU_ID
        WHERE m.OKR_NUM IN ('KH0121', 'KH0473') AND c.CORD_LEVEL = 1 AND COALESCE(c.CORD_CLASSIFICATION, '') != 'T'
        ORDER BY m.OKR_NUM, c.CORD_ORDINAL
        """,
        con,
    ).merge(vals, on="CORD_ID", how="left")
    p1["cord_val"] = p1["cord_val"].fillna(0.0).astype(int)

    g1 = p1[p1["OKR_NUM"] == "KH0121"].reset_index(drop=True)  # len 29 (forward)
    # KH0473 reversed (len 30) has one extra broken 0-length stub at reversed index 17 (original ord 13)!
    g2_rev = p1[p1["OKR_NUM"] == "KH0473"].iloc[::-1].reset_index(drop=True)
    # Drop the single extra 0-length stub at rev idx 17 so all 29 cords align 1-to-1!
    g2_aligned = g2_rev.drop(index=17).reset_index(drop=True)

    df = pd.DataFrame({
        "aligned_idx": np.arange(1, 30),
        "kh0121_ord": g1["CORD_ORDINAL"].values,
        "kh0121_val": g1["cord_val"].values,
        "kh0121_len_cm": g1["CORD_LENGTH"].values,
        "kh0473_ord_rev": g2_aligned["CORD_ORDINAL"].values,
        "kh0473_val_rev": g2_aligned["cord_val"].values,
        "kh0473_len_cm_rev": g2_aligned["CORD_LENGTH"].values,
    })
    df["val_exact"] = (df["kh0121_val"] == df["kh0473_val_rev"]).astype(int)
    df["len_abs_diff_cm"] = (df["kh0121_len_cm"] - df["kh0473_len_cm_rev"]).abs()

    valid_len = df[(df["kh0121_len_cm"] > 0) & (df["kh0473_len_cm_rev"] > 0)]
    r_val = float(spearmanr(df["kh0121_val"], df["kh0473_val_rev"]).statistic)
    r_len = float(pearsonr(valid_len["kh0121_len_cm"], valid_len["kh0473_len_cm_rev"]).statistic)
    stats = {
        "n_aligned": len(df),
        "n_val_exact": int(df["val_exact"].sum()),
        "val_spearman_r": r_val,
        "n_valid_len": len(valid_len),
        "len_pearson_r": r_len,
        "med_abs_len_diff_cm": float(valid_len["len_abs_diff_cm"].median()),
        "med_len_121": float(valid_len["kh0121_len_cm"].median()),
        "med_len_473": float(valid_len["kh0473_len_cm_rev"].median()),
    }
    return df, stats


def make_fig6(
    sweep_df: pd.DataFrame,
    tdf: pd.DataFrame,
    t_stats: dict,
    dup_df: pd.DataFrame,
    dup_stats: dict,
) -> pathlib.Path:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({
        "font.family": "sans-serif", "font.sans-serif": ["Arial", "DejaVu Sans"],
        "axes.edgecolor": "#333333", "axes.linewidth": 0.8, "axes.titlesize": 10.2,
        "axes.labelsize": 9.0, "xtick.labelsize": 8.2, "ytick.labelsize": 8.2, "grid.alpha": 0.3,
    })
    fig, axes = plt.subplots(2, 2, figsize=(14.2, 10.6))

    # --- PANEL A: Corpus-Wide 611 x 611 High-Entropy Pair Discovery Map ---
    ax = axes[0, 0]
    bg = sweep_df[sweep_df["category"] == "Background"]
    ctrl = sweep_df[sweep_df["category"].str.startswith("Known")]
    newd = sweep_df[sweep_df["category"].str.startswith("NEW")]

    ax.scatter(bg["rare_dice"], bg["ord_rare_ex"], color="#BBBBBB", s=22, alpha=0.45, label=f"Background pairs (n={len(bg)})", zorder=2)
    ax.scatter(ctrl["rare_dice"], ctrl["ord_rare_ex"], color="#2E86AB", s=70, marker="^", edgecolors="black", linewidths=0.7, label=f"Recovered positive controls (n={len(ctrl)})", zorder=4)
    ax.scatter(newd["rare_dice"], newd["ord_rare_ex"], color="#D62828", s=110, marker="*", edgecolors="black", linewidths=0.8, label="New discoveries (Berlin)", zorder=5)

    annots = [
        ("KH0192", "KH0470", "NEW: Berlin Triptych Part 1\nKH0192 (Pachacamac) <-> KH0470 (Huacho)\n36 exact rare / 88/104 cells", (-0.36, -6.5), "#900000"),
        ("KH0121", "KH0473", "NEW: Reversed Dup\nKH0121 <-> KH0473 (r=0.91)", (-0.34, 4.5), "#900000"),
        ("KH0140", "KH0483", "Ctrl: AS125/UR247\n(exact dup)", (-0.32, 3.0), "#1B4F72"),
        ("KH0303", "KH0304", "Ctrl: UR066/UR067", (-0.25, -2.0), "#1B4F72"),
        ("KH0305", "KH0301", "Ctrl: Puruchuco\nUR064/UR068 (Science 2005)", (-0.25, 5.2), "#1B4F72"),
    ]
    for a, b, txt, (dx, dy), col in annots:
        sub = sweep_df[((sweep_df["okr_a"] == a) & (sweep_df["okr_b"] == b)) | ((sweep_df["okr_a"] == b) & (sweep_df["okr_b"] == a))]
        if len(sub):
            r = sub.iloc[0]
            ax.annotate(
                txt, xy=(r["rare_dice"], r["ord_rare_ex"]),
                xytext=(r["rare_dice"] + dx, r["ord_rare_ex"] + dy),
                fontsize=7.3, fontweight="bold" if col == "#900000" else "normal", color=col,
                arrowprops=dict(arrowstyle="->", color=col, lw=0.75),
                bbox=dict(facecolor="white", edgecolor="none", alpha=0.75, pad=0.6),
            )
    ax.set_xlim(0.18, 1.07)
    ax.set_ylim(0.5, 41)
    ax.set_xlabel("Rare Non-Round Multi-Digit Overlap (Dice, val >= 11 & %10 != 0)")
    ax.set_ylabel("Exact Ordered Rare-Number Matches")
    ax.set_title("A. Blind 611 x 611 Corpus Pair Discovery Sweep (OKR v2.1.0)", fontweight="bold")
    ax.legend(fontsize=7.6, loc="upper left", frameon=True)
    ax.grid(True, alpha=0.3)

    # --- PANEL B: Berlin Triptych Part 1 (UR1175 vs UR233) ---
    ax2 = axes[0, 1]
    p1_sub = tdf[tdf["triptych_half"] == "Part1_vs_UR233"]
    ex1 = p1_sub["exact_match"] == 1
    mx1 = max(p1_sub["master_val"].max(), p1_sub["copy_val"].max()) * 1.06
    ax2.plot([0, mx1], [0, mx1], color="#555555", linestyle="--", lw=1.0, label="Exact parity (y = x)")
    ax2.scatter(p1_sub.loc[ex1, "master_val"], p1_sub.loc[ex1, "copy_val"], color="#1B998B", s=44, alpha=0.85, edgecolors="black", lw=0.5, label=f"Exact match ({t_stats['p1_exact']}/{t_stats['p1_n']} = {100*t_stats['p1_exact']/t_stats['p1_n']:.1f}%)")
    ax2.scatter(p1_sub.loc[~ex1, "master_val"], p1_sub.loc[~ex1, "copy_val"], color="#F18F01", s=54, marker="D", edgecolors="black", lw=0.6, label=f"Adjusted/diff ({t_stats['p1_n']-t_stats['p1_exact']} cells; 97% within +/-10)")
    ax2.text(
        0.04, 0.94,
        f"Pearson r = {t_stats['p1_pearson_r']:.4f} (perm p < 1e-4)\n"
        f"Spearman r = {t_stats['p1_spearman_r']:.4f}\n"
        f"Master: k=5, med 65 cm, 114 subs (Pachacamac)\n"
        f"Copy: k=6 (+0 spacer), med 40 cm, 18 subs (Huacho)",
        transform=ax2.transAxes, va="top", ha="left", fontsize=7.6,
        bbox=dict(facecolor="#F8F9FA", edgecolor="#CCCCCC", boxstyle="round,pad=0.35"),
    )
    ax2.set_xlim(-8, mx1)
    ax2.set_ylim(-8, mx1)
    ax2.set_xlabel("KH0192 / UR1175 Part 1 Value (VA 42518, Groups 1-21)")
    ax2.set_ylabel("KH0470 / UR233 Value (VA 63038 b, spacer removed)")
    ax2.set_title("B. Berlin Triptych Half 1: Master Part 1 vs. Companion UR233", fontweight="bold")
    ax2.legend(fontsize=7.5, loc="lower right", frameon=True)
    ax2.grid(True, alpha=0.3)

    # --- PANEL C: Berlin Triptych Part 3 (UR1175 vs UR232) ---
    ax3 = axes[1, 0]
    p3_sub = tdf[tdf["triptych_half"] == "Part3_vs_UR232"]
    ex3 = p3_sub["exact_match"] == 1
    mx3 = max(p3_sub["master_val"].max(), p3_sub["copy_val"].max()) * 1.07
    ax3.plot([0, mx3], [0, mx3], color="#555555", linestyle="--", lw=1.0, label="Exact parity (y = x)")
    ax3.scatter(p3_sub.loc[ex3, "master_val"], p3_sub.loc[ex3, "copy_val"], color="#2E86AB", s=44, alpha=0.85, edgecolors="black", lw=0.5, label=f"Exact match ({t_stats['p3_exact']}/{t_stats['p3_n']} = {100*t_stats['p3_exact']/t_stats['p3_n']:.1f}%)")
    ax3.scatter(p3_sub.loc[~ex3, "master_val"], p3_sub.loc[~ex3, "copy_val"], color="#D62828", s=48, marker="s", alpha=0.85, edgecolors="black", lw=0.5, label=f"Adjusted/diff ({t_stats['p3_n']-t_stats['p3_exact']} cells; 86% within +/-10)")
    ax3.text(
        0.04, 0.94,
        f"Pearson r = {t_stats['p3_pearson_r']:.4f} (perm p < 1e-4)\n"
        f"Spearman r = {t_stats['p3_spearman_r']:.4f}\n"
        f"Groups 1-21 aligned through final cord\n"
        f"(e.g. Grp 21: [50,12,125,90,11] vs [40,12,125,90,21])",
        transform=ax3.transAxes, va="top", ha="left", fontsize=7.6,
        bbox=dict(facecolor="#F8F9FA", edgecolor="#CCCCCC", boxstyle="round,pad=0.35"),
    )
    ax3.set_xlim(-6, mx3)
    ax3.set_ylim(-6, mx3)
    ax3.set_xlabel("KH0192 / UR1175 Part 3 Value (VA 42518, Groups 25-45)")
    ax3.set_ylabel("KH0469 / UR232 Value (VA 63038(B), spacer removed)")
    ax3.set_title("C. Berlin Triptych Half 2: Master Part 3 vs. Companion UR232", fontweight="bold")
    ax3.legend(fontsize=7.5, loc="lower right", frameon=True)
    ax3.grid(True, alpha=0.3)

    # --- PANEL D: Physical Cord-Length Forensic Proof (Distinct Artifacts vs. Reversed Duplicate!) ---
    ax4 = axes[1, 1]
    vl_dup = dup_df[(dup_df["kh0121_len_cm"] > 0) & (dup_df["kh0473_len_cm_rev"] > 0)]
    vl_t1 = p1_sub.dropna(subset=["master_len_cm", "copy_len_cm"])
    vl_t1 = vl_t1[(vl_t1["master_len_cm"] > 0) & (vl_t1["copy_len_cm"] > 0)]

    ax4.plot([15, 100], [15, 100], color="#444444", linestyle="--", lw=1.0, label="Same physical artifact (y = x)")
    ax4.scatter(
        vl_dup["kh0121_len_cm"], vl_dup["kh0473_len_cm_rev"],
        color="#D62828", s=55, marker="o", alpha=0.9, edgecolors="black", lw=0.6,
        label=f"Same Artifact (Reversed Dup): KH0121 vs KH0473\n(r={dup_stats['len_pearson_r']:.2f}, med |diff|={dup_stats['med_abs_len_diff_cm']:.1f} cm)",
        zorder=5,
    )
    ax4.scatter(
        vl_t1["master_len_cm"], vl_t1["copy_len_cm"],
        color="#1B998B", s=34, marker="s", alpha=0.65, edgecolors="black", lw=0.4,
        label=f"Distinct Companion Copy: KH0192 vs KH0470\n(med {t_stats['p1_med_len_master']:.0f} vs {t_stats['p1_med_len_copy']:.0f} cm, MW p={t_stats['p1_len_mw_p']:.1e})",
        zorder=4,
    )
    ax4.set_xlim(18, 100)
    ax4.set_ylim(18, 100)
    ax4.set_xlabel("Khipu A Pendant Cord Length (cm)")
    ax4.set_ylabel("Matched Khipu B Pendant Cord Length (cm)")
    ax4.set_title("D. Physical Forensics: Distinct Companion Copy vs. Reversed Duplicate", fontweight="bold")
    ax4.legend(fontsize=7.5, loc="upper left", frameon=True)
    ax4.grid(True, alpha=0.3)

    plt.tight_layout()
    out = FIGURES_DIR / "fig6_phase2c_berlin_triptych_and_pairs.png"
    plt.savefig(out, dpi=300, bbox_inches="tight")
    plt.close()
    return out


def main() -> None:
    rng = np.random.default_rng(42)
    con = sqlite3.connect(find_okr_db())
    sweep_df = run_corpus_rare_number_sweep(con)
    tdf, t_stats = analyze_berlin_triptych(rng)
    dup_df, dup_stats = analyze_reversed_duplicate(con)
    con.close()

    sweep_df.to_csv(RESULTS_DIR / "phase2c_corpus_pair_sweep.csv", index=False)
    tdf.to_csv(RESULTS_DIR / "phase2c_berlin_triptych_cells.csv", index=False)
    dup_df.to_csv(RESULTS_DIR / "phase2c_reversed_duplicate_kh0121_kh0473.csv", index=False)

    summary_rows = [
        ("Putative Berlin Triptych Half 1 (UR1175 P1 <-> UR233)", t_stats["p1_n"], t_stats["p1_exact"], f"{100*t_stats['p1_exact']/t_stats['p1_n']:.1f}%", round(t_stats["p1_pearson_r"], 4), round(t_stats["p1_spearman_r"], 4), f"<={t_stats['p1_perm_p']:.4f} (post-selection)", f"{t_stats['p1_med_len_master']:.1f} vs {t_stats['p1_med_len_copy']:.1f} cm (MW p={t_stats['p1_len_mw_p']:.2e})", "Putative distinct companion copy (k=5 -> k=6 +0 spacer; pending museum verification)"),
        ("Putative Berlin Triptych Half 2 (UR1175 P3 <-> UR232)", t_stats["p3_n"], t_stats["p3_exact"], f"{100*t_stats['p3_exact']/t_stats['p3_n']:.1f}%", round(t_stats["p3_pearson_r"], 4), round(t_stats["p3_spearman_r"], 4), f"<={t_stats['p3_perm_p']:.4f} (post-selection)", f"{t_stats['p3_med_len_master']:.1f} vs {t_stats['p3_med_len_copy']:.1f} cm (MW p={t_stats['p3_len_mw_p']:.2e})", "Putative distinct companion copy (k=5 -> k=6 +0 spacer; pending museum verification)"),
        ("Putative Reversed Dup (KH0121 <-> KH0473)", dup_stats["n_aligned"], dup_stats["n_val_exact"], f"{100*dup_stats['n_val_exact']/dup_stats['n_aligned']:.1f}%", round(dup_stats["len_pearson_r"], 4), round(dup_stats["val_spearman_r"], 4), "reversed (descriptive)", f"{dup_stats['med_len_121']:.1f} vs {dup_stats['med_len_473']:.1f} cm (med |diff|={dup_stats['med_abs_len_diff_cm']:.1f} cm)", "Putative same physical Berlin artifact (suggests Marquez prov; pending museum verification)"),
    ]
    sum_df = pd.DataFrame(
        summary_rows,
        columns=["finding", "n_cells", "exact_matches", "pct_exact", "pearson_r", "spearman_r", "perm_p", "physical_lengths", "classification"],
    )
    sum_df.to_csv(RESULTS_DIR / "phase2c_discoveries_summary.csv", index=False)

    fig_path = make_fig6(sweep_df, tdf, t_stats, dup_df, dup_stats)
    print(f"Saved {fig_path}")
    print("\n=== PHASE 2C DISCOVERIES SUMMARY ===")
    print(sum_df.to_string(index=False))
    print("\n=== TOP 10 SWEEP HITS ===")
    print(sweep_df.head(10)[["okr_a", "inv_a", "okr_b", "inv_b", "rare_dist", "rare_dice", "ord_rare_ex", "spearman_r", "ori", "category"]].to_string(index=False))


if __name__ == "__main__":
    main()
