"""
step2_taxonomy_outliers.py — Pre-registered Step 2:
Full structural anatomy, deterministic family taxonomy (F1-F4 + OUT),
unforced HDBSCAN noise accounting, and individual outlier dossier for the 86 doubly-tabular khipus.

Uses bug-free subsidiary counts traced directly from OKR cord table (never KFG khipu_summary.csv)
and cross-checks exact x100 group totals via KFG_DB_Excel workbooks.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.cluster import HDBSCAN
from sklearn.preprocessing import StandardScaler

from src.data_loader import (
    RESULTS_DIR,
    extract_modal_table,
    load_cohorts,
    load_kfg_excel_pendants,
    load_okr_pendants,
)


def compute_x100_groups(mat_val: np.ndarray, kfg_df: pd.DataFrame | None = None) -> tuple[int, int]:
    """Count groups (>=2 non-zero cords) summing to exact multiples of 100 (loose and strict)."""
    # Prefer KFG Excel clusters when available for exact agreement with v2 KFG analysis
    if kfg_df is not None and len(kfg_df) > 0:
        valid = kfg_df[kfg_df["kfg_cluster"] >= 0]
        loose = strict = 0
        for _, grp in valid.groupby("kfg_cluster"):
            v = grp["kfg_val"].values.astype(int)
            s = int(v.sum())
            if s >= 100 and s % 100 == 0 and int(np.sum(v > 0)) >= 2:
                loose += 1
                if int(np.sum((v % 10) != 0)) >= 1:
                    strict += 1
        return loose, strict

    loose = strict = 0
    for row in mat_val:
        v = np.round(row).astype(int)
        s = int(v.sum())
        if s >= 100 and s % 100 == 0 and int(np.sum(v > 0)) >= 2:
            loose += 1
            if int(np.sum((v % 10) != 0)) >= 1:
                strict += 1
    return loose, strict


def subsidiary_slot_concentration(mat_sub: np.ndarray) -> tuple[int, float, float]:
    """
    Test whether subsidiaries attach preferentially to specific column slots 1..k.
    Returns (total_subs_in_table, normalized_entropy [0..1], chi2_p_value).
    """
    col_sums = mat_sub.sum(axis=0).astype(float)
    tot = int(col_sums.sum())
    k = len(col_sums)
    if tot == 0 or k < 2:
        return tot, np.nan, np.nan
    p = col_sums / tot
    nz = p > 0
    ent = -float(np.sum(p[nz] * np.log2(p[nz]))) / np.log2(k)
    if tot >= 5:
        _, pval = stats.chisquare(col_sums)
        pval = float(pval)
    else:
        pval = np.nan
    return tot, round(ent, 4), pval


def binary_purity(mat: np.ndarray, pos_tok: str, neg_tok: str) -> tuple[float, float, float]:
    """
    Compute:
      - known_frac: fraction of cells with pos_tok or neg_tok
      - row_purity: mean within-row majority share across rows with >=2 known cords
      - overall_pos_share: share of pos_tok among known cords
    """
    is_pos = mat == pos_tok
    is_neg = mat == neg_tok
    known = is_pos | is_neg
    n_known = int(known.sum())
    if n_known == 0:
        return 0.0, np.nan, np.nan
    known_frac = n_known / mat.size
    overall_pos = float(is_pos.sum() / n_known)

    row_purities = []
    for r in range(mat.shape[0]):
        rk = known[r]
        if rk.sum() >= 2:
            p = is_pos[r][rk].mean()
            row_purities.append(max(p, 1.0 - p))
    row_pur = float(np.mean(row_purities)) if row_purities else np.nan
    return round(known_frac, 4), (round(row_pur, 4) if np.isfinite(row_pur) else np.nan), round(overall_pos, 4)


def assign_deterministic_family(row: pd.Series) -> tuple[str, str]:
    """
    Assign pre-registered deterministic families F1-F4 or OUT (outlier) with explicit reason.
    Strictly follows PREREGISTRATION.md Section 2.
    """
    k = int(row["k"])
    G = int(row["G"])
    clean = float(row["clean_dec_frac"])
    h5c_sig = bool(row["h5c_sig"])
    x100_strict = int(row["n_groups_x100_strict"])

    reasons = []
    if k >= 10:
        reasons.append(f"extreme_arity(k={k})")
    if G == 2:
        reasons.append("minimal_rows(G=2)")
    if row["table_coverage"] < 0.25:
        reasons.append(f"low_table_coverage({row['table_coverage']:.2f})")

    # Extreme arity (k >= 10 except F3 self-totaling up to k=12) or G=2 -> OUT
    is_self_tot = h5c_sig or (x100_strict >= 2)
    if (3 <= k <= 12) and G >= 3 and clean >= 0.90 and is_self_tot:
        return "F3_self_totaling", ""
    if k >= 10 or G < 3:
        return "OUT", "; ".join(reasons) if reasons else f"outside_F1_F4(G={G},k={k})"
    if clean < 0.90:
        return "F4_nonclean_hybrid", ""
    if k in (3, 4):
        return "F1_compact_parallel", ""
    if 5 <= k <= 9:
        return "F2_medium_multi_attr", ""
    return "OUT", "; ".join(reasons) if reasons else "unclassified"


def main():
    print("Loading OKR v2.1.0 pendants, subsidiaries, KFG Excel values, and cohorts...")
    pendants, khipu_sub_totals = load_okr_pendants()
    cohorts = load_cohorts()
    dbl = cohorts[cohorts["doubly_tabular"]].copy().reset_index(drop=True)
    print(f"Doubly-tabular cohort size: {len(dbl)} (expected 86)")
    assert len(dbl) == 86, f"Expected 86 doubly-tabular khipus, got {len(dbl)}"

    pend_by_kid = {int(kid): g for kid, g in pendants.groupby("KHIPU_ID", sort=False)}
    sig_rows = []

    for r in dbl.itertuples(index=False):
        kid = int(r.KHIPU_ID)
        g = pend_by_kid[kid]
        tab = extract_modal_table(g)
        assert tab is not None, f"Missing modal table for {r.OKR_NUM}"

        G, k = tab["G"], tab["k"]
        n_pend = len(g)
        cov = (G * k) / n_pend
        aspect = G / k

        log_mat = np.log1p(tab["mat_val"])
        col_means = log_mat.mean(axis=0)
        mono_rho = float(stats.spearmanr(np.arange(k), col_means).statistic) if np.std(col_means) > 1e-12 else 0.0
        col_zero_rates = (tab["mat_val"] == 0).mean(axis=0)
        col_zero_var = float(np.var(col_zero_rates))

        # Check KFG Excel for exact x100 totals
        lf = r.legal_filename if isinstance(r.legal_filename, str) else ""
        kfg_df = load_kfg_excel_pendants(lf) if lf else None
        x100_loose, x100_strict = compute_x100_groups(tab["mat_val"], kfg_df)

        h5c_sig = bool(
            (pd.notna(r.h5c_p) and r.h5c_p < 0.05) or (pd.notna(r.h5c_any_p) and r.h5c_any_p < 0.05)
        )

        # Bug-free subsidiaries from OKR
        n_sub_total = int(khipu_sub_totals.get(kid, 0))
        sub_ratio = n_sub_total / max(n_pend, 1)
        sub_in_tab, sub_ent, sub_chi2_p = subsidiary_slot_concentration(tab["mat_sub"])

        # Binary axes
        rv_known, rv_row_pur, rv_r_share = binary_purity(tab["mat_att"], "R", "V")
        tw_known, tw_row_pur, tw_s_share = binary_purity(tab["mat_twist"], "S", "Z")
        ks_tot = int(tab["mat_ks"].sum())
        kz_tot = int(tab["mat_kz"].sum())
        knot_s_share = ks_tot / (ks_tot + kz_tot) if (ks_tot + kz_tot) > 0 else np.nan

        sig_rows.append({
            "KHIPU_ID": kid,
            "OKR_NUM": r.OKR_NUM,
            "kfg_name": r.kfg_name if pd.notna(r.kfg_name) else "",
            "cataloger": r.cataloger,
            "provenance": r.provenance_clean,
            "museum": r.museum_clean,
            "keep_dedup": bool(r.keep_dedup),
            "n_pendants": n_pend,
            "G": G,
            "k": k,
            "table_coverage": round(cov, 4),
            "aspect_ratio": round(aspect, 4),
            "clean_dec_frac": round(float(r.clean_dec_frac), 4),
            "h3_col_match": round(float(r.h3_col_match), 4),
            "h3_val_r2": round(float(r.h3_val_r2), 4),
            "col_mono_rho": round(mono_rho, 4),
            "abs_col_mono_rho": round(abs(mono_rho), 4),
            "col_zero_var": round(col_zero_var, 4),
            "n_groups_x100": x100_loose,
            "n_groups_x100_strict": x100_strict,
            "h5c_hits": r.h5c_hits,
            "h5c_p": r.h5c_p,
            "h5c_any_hits": r.h5c_any_hits,
            "h5c_any_p": r.h5c_any_p,
            "h5c_sig": h5c_sig,
            "n_sub_total": n_sub_total,
            "sub_ratio": round(sub_ratio, 4),
            "sub_in_table": sub_in_tab,
            "sub_col_entropy": sub_ent,
            "sub_col_chi2_p": round(sub_chi2_p, 5) if np.isfinite(sub_chi2_p) else np.nan,
            "rv_known_frac": rv_known,
            "rv_row_purity": rv_row_pur,
            "rv_r_share": rv_r_share,
            "twist_known_frac": tw_known,
            "twist_s_share": tw_s_share,
            "knot_s_share": round(knot_s_share, 4) if np.isfinite(knot_s_share) else np.nan,
        })

    sig_df = pd.DataFrame(sig_rows)

    # 1. Deterministic families + explicit OUT
    fams = [assign_deterministic_family(row) for _, row in sig_df.iterrows()]
    sig_df["det_family"] = [f[0] for f in fams]
    sig_df["det_outlier_reason"] = [f[1] for f in fams]
    sig_df["is_det_outlier"] = sig_df["det_family"] == "OUT"

    # 2. Secondary unsupervised check: HDBSCAN with unforced noise (-1)
    hdb_feats = [
        "log_G", "log_k", "table_coverage", "abs_col_mono_rho",
        "col_zero_var", "sub_ratio", "h3_col_match", "h3_val_r2", "clean_dec_frac",
    ]
    X_raw = sig_df.assign(
        log_G=np.log1p(sig_df["G"]),
        log_k=np.log1p(sig_df["k"]),
    )[hdb_feats].fillna(0.0).values
    X_std = StandardScaler().fit_transform(X_raw)
    hdb = HDBSCAN(min_cluster_size=5, min_samples=3)
    sig_df["hdbscan_cluster"] = hdb.fit_predict(X_std)
    sig_df["is_hdbscan_noise"] = sig_df["hdbscan_cluster"] == -1

    # Combined outlier indicator for dossier
    sig_df["is_any_outlier"] = sig_df["is_det_outlier"] | sig_df["is_hdbscan_noise"]

    # Save full signature table
    sig_df.to_csv(RESULTS_DIR / "step2_signatures_86.csv", index=False)

    # Family summary table
    fam_sum = sig_df.groupby("det_family").agg(
        n=("KHIPU_ID", "count"),
        pct=("KHIPU_ID", lambda s: round(100.0 * len(s) / len(sig_df), 1)),
        med_G=("G", "median"),
        med_k=("k", "median"),
        med_cov=("table_coverage", "median"),
        med_col_match=("h3_col_match", "median"),
        med_val_r2=("h3_val_r2", "median"),
        med_sub_ratio=("sub_ratio", "median"),
        n_hdb_noise=("is_hdbscan_noise", "sum"),
    ).reset_index()
    fam_sum.to_csv(RESULTS_DIR / "step2_family_summary.csv", index=False)

    # Subsidiary slot-concentration summary (BH-FDR among khipus with >=5 table subsidiaries)
    sub_testable = sig_df[sig_df["sub_col_chi2_p"].notna()].copy()
    if len(sub_testable) > 0:
        from scipy.stats import false_discovery_control
        sub_testable["sub_col_q"] = false_discovery_control(sub_testable["sub_col_chi2_p"].values, method="bh")
        sig_df = sig_df.merge(sub_testable[["KHIPU_ID", "sub_col_q"]], on="KHIPU_ID", how="left")
        sig_df.to_csv(RESULTS_DIR / "step2_signatures_86.csv", index=False)

    # Outlier dossier (all deterministic OUT + all HDBSCAN noise)
    dossier_cols = [
        "OKR_NUM", "kfg_name", "provenance", "cataloger", "G", "k", "n_pendants",
        "table_coverage", "clean_dec_frac", "h3_col_match", "h3_val_r2",
        "n_sub_total", "sub_ratio", "h5c_sig", "n_groups_x100_strict",
        "det_family", "det_outlier_reason", "hdbscan_cluster",
    ]
    dossier = sig_df[sig_df["is_any_outlier"]].sort_values(
        ["is_det_outlier", "k", "G"], ascending=[False, False, False]
    )[dossier_cols]
    dossier.to_csv(RESULTS_DIR / "step2_outliers_dossier.csv", index=False)

    n_det_out = int(sig_df["is_det_outlier"].sum())
    n_hdb_noise = int(sig_df["is_hdbscan_noise"].sum())
    n_both_out = int((sig_df["is_det_outlier"] & sig_df["is_hdbscan_noise"]).sum())
    n_hdb_clusters = len(set(sig_df["hdbscan_cluster"]) - {-1})

    print("\n=== STEP 2: DETERMINISTIC FAMILY SUMMARY (n = 86) ===")
    print(fam_sum.to_string(index=False))

    print("\n=== STEP 2: EXPLICIT OUTLIER ACCOUNTING ===")
    print(f"  Deterministic OUT (outside F1-F4): {n_det_out} / 86 ({100*n_det_out/86:.1f}%)")
    print(f"  HDBSCAN unforced noise (-1):       {n_hdb_noise} / 86 ({100*n_hdb_noise/86:.1f}%)  [HDBSCAN found {n_hdb_clusters} clusters]")
    print(f"  Overlap (both OUT & HDBSCAN -1):   {n_both_out} / 86 ({100*n_both_out/86:.1f}%)")

    if len(sub_testable) > 0:
        n_sub_sig = int((sub_testable["sub_col_q"] < 0.05).sum())
        print(f"\n=== STEP 2: SUBSIDIARY COLUMN CONCENTRATION (bug-free OKR) ===")
        print(f"  Testable (>=5 subs in modal table): {len(sub_testable)} / 86")
        print(f"  Non-uniform slot attachment (BH q<0.05): {n_sub_sig} / {len(sub_testable)} ({100*n_sub_sig/len(sub_testable):.1f}%)")
        print(f"  Median normalized slot entropy: {sub_testable['sub_col_entropy'].median():.3f} (1.0 = uniform)")

    print("\n=== STEP 2: DETERMINISTIC OUTLIERS DOSSIER (det_family == 'OUT') ===")
    det_out_df = dossier[dossier["det_family"] == "OUT"]
    print(det_out_df[["OKR_NUM", "kfg_name", "provenance", "G", "k", "n_pendants", "table_coverage", "clean_dec_frac", "n_sub_total", "det_outlier_reason", "hdbscan_cluster"]].to_string(index=False))


if __name__ == "__main__":
    main()
