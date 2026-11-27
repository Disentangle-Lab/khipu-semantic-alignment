"""
step1_2_hardening.py — Statistical & archaeological hardening block (#2, #4, #5):
  - #2 (Step 1 Audit): Exact FULL_COLOR sequences, museum numbers, tomb/bundle metadata
    for every (site, k) cell in primary_86, plus plain-cotton vs complex/dyed color check.
  - #4 (Step 2 Subsidiaries): Replace asymptotic chi-square with exact intra-row
    permutation test (N_PERM = 5,000) on BOTH raw subsidiary counts AND binary
    has-subsidiary indicators (immune to within-cord overdispersion & small counts).
  - #5 (Step 2 Outliers): Exact Fisher tests for enrichment of deterministic OUT (21/86)
    and HDBSCAN noise (45/86) in (a) official calendar khipus, (b) S/Z spatial-block
    calendar candidates, and (c) KFG Part B newly detected paired khipus.
"""

from __future__ import annotations

import sqlite3
from collections import Counter

import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import false_discovery_control

from src.data_loader import (
    RESULTS_DIR,
    extract_modal_table,
    find_okr_db,
    load_cohorts,
    load_okr_pendants,
)
from src.step1_archive_invariance import modal_color_profile

N_PERM = 5000
SEED = 2026

# Plain undyed cotton / neutral brown/white codes in Ascher color coding
PLAIN_COTTON_CODES = {"W", "AB", "MB", "LB", "DB", "YB", "KB", "B", "GB", "EB", "FB", "BB"}

# Official KFG S/Z calendar candidates (from andean-khipus-kfg-final.ipynb) + UR006/KH0242 Calendar Khipu
OFFICIAL_CALENDAR_KFG10 = {
    "KH0300", "KH0311", "KH0033", "KH0063", "KH0097",
    "KH0118", "KH0349", "KH0405", "KH0433", "KH0436",
}
OFFICIAL_CALENDAR_PLUS_UR006 = OFFICIAL_CALENDAR_KFG10 | {"KH0242"}

# Exact 9 new paired-khipu candidates from kfg-khipu-analyses.ipynb Part B (kfg_name / okr_name)
PARTB_NEW_PAIR_NAMES = {
    "AS086", "MM1086",
    "AS125", "UR247",
    "AS067/MA029", "AS067", "KH0080",
    "UR066", "UR067",
    "UR064", "UR068",
    "UR255", "UR267A",
    "HP055", "UR1084",
    "UR256", "UR267B",
    "UR063", "UR073",
}


def is_complex_color(code: str) -> bool:
    if not code or code == "__MISSING__":
        return False
    if any(sep in code for sep in (":", "-", "/", "*", "@", "%")):
        return True
    return code not in PLAIN_COTTON_CODES


def run_part_a_step1_audit(pendants: pd.DataFrame, cohorts: pd.DataFrame) -> pd.DataFrame:
    """Inspect every khipu in primary_86 (site, k) cells: museum_num, notes, exact colors, numeric profile."""
    con = sqlite3.connect(str(find_okr_db()))
    km = pd.read_sql(
        "SELECT KHIPU_ID, MUSEUM_NAME, MUSEUM_NUM, DATE_DISCOVERED, DISCOVERED_BY, ARCHIVE_NUM, NOTES, NICKNAME FROM khipu_main",
        con,
    )
    con.close()

    dbl = cohorts[cohorts["doubly_tabular"]].merge(km, on="KHIPU_ID", how="left", suffixes=("", "_km"))
    pend_by_kid = {int(k): g for k, g in pendants.groupby("KHIPU_ID", sort=False)}

    # Identify all (site, k) cells with n >= 2 (excluding Unknown)
    known = dbl[dbl["provenance_clean"] != "Unknown"].copy()
    cell_counts = known.groupby(["provenance_clean", "h3_m"]).size()
    valid_cells = set(cell_counts[cell_counts >= 2].index.tolist())

    rows = []
    for r in known.itertuples(index=False):
        key = (r.provenance_clean, int(r.h3_m))
        if key not in valid_cells:
            continue
        kid = int(r.KHIPU_ID)
        g = pend_by_kid[kid]
        tab = extract_modal_table(g)
        mc = modal_color_profile(tab["mat_col"])
        col_means = np.log1p(tab["mat_val"]).mean(axis=0)
        ranks = stats.rankdata(col_means).round(1).tolist()
        n_complex = sum(is_complex_color(c) for c in mc)
        notes_short = str(r.NOTES or "").replace("\n", " ").replace("\r", " ")[:140]
        rows.append({
            "site": r.provenance_clean,
            "k": int(r.h3_m),
            "G": int(r.h3_G),
            "OKR_NUM": r.OKR_NUM,
            "kfg_name": r.kfg_name if pd.notna(r.kfg_name) else "",
            "cataloger": r.cataloger,
            "keep_dedup": bool(r.keep_dedup),
            "museum": str(r.museum_clean)[:45],
            "museum_num": str(r.MUSEUM_NUM or ""),
            "archive_num": r.ARCHIVE_NUM,
            "discovered_by": str(r.DISCOVERED_BY or ""),
            "modal_colors": " | ".join(mc),
            "n_complex_cols": n_complex,
            "frac_complex_cols": round(n_complex / len(mc), 3),
            "col_mean_log1p": ", ".join(f"{v:.2f}" for v in col_means),
            "col_ranks": str(ranks),
            "notes_snippet": notes_short,
        })

    audit_df = pd.DataFrame(rows).sort_values(["site", "k", "OKR_NUM"]).reset_index(drop=True)
    audit_df.to_csv(RESULTS_DIR / "step1_cells_detailed_audit.csv", index=False)

    print("=== PART A (#2): EXACT COLORS & MUSEUM/TOMB METADATA BY (SITE, k) CELL ===")
    for (site, k), grp in audit_df.groupby(["site", "k"]):
        print(f"\n--- {site} | k={k} (n={len(grp)}) ---")
        # Compute slot-wise modal color sharing across khipus in this cell
        seqs = [s.split(" | ") for s in grp["modal_colors"]]
        all_codes = [c for seq in seqs for c in seq if c != "__MISSING__"]
        top_c = Counter(all_codes).most_common(6)
        cx_rate = np.mean([is_complex_color(c) for c in all_codes]) if all_codes else 0.0
        print(f"  Top colors in cell: {top_c} | Complex/dyed share: {100*cx_rate:.1f}%")
        for r in grp.itertuples(index=False):
            print(
                f"  {r.OKR_NUM} ({r.kfg_name:7s}, {r.cataloger}, dedup={int(r.keep_dedup)}) "
                f"mus_num={r.museum_num!r:20s} | cols=[{r.modal_colors}] | ranks={r.col_ranks}"
            )
            if r.notes_snippet:
                print(f"      notes: {r.notes_snippet}")
    return audit_df


def run_part_b_subsidiary_perm(pendants: pd.DataFrame, sig86: pd.DataFrame) -> pd.DataFrame:
    """
    Exact intra-row permutation test (N_PERM = 5,000) for subsidiary column concentration:
      1) Raw counts: mat_sub (G x k), permuting slots within each row g independently
         (preserves row totals AND within-cord subsidiary clumping/overdispersion!).
         Statistic: normalized entropy H_norm (lower = more concentrated) and
         slot-sum variance V = sum((col_sum - mean)^2) (higher = more concentrated).
      2) Binary presence: mat_bin = (mat_sub > 0).astype(int), same intra-row permutation
         (completely immune to a single pendant carrying 20 subsidiaries!).
    """
    rng = np.random.default_rng(SEED)
    pend_by_kid = {int(k): g for k, g in pendants.groupby("KHIPU_ID", sort=False)}
    rows = []

    for r in sig86.itertuples(index=False):
        kid = int(r.KHIPU_ID)
        tab = extract_modal_table(pend_by_kid[kid])
        mat_sub = tab["mat_sub"].astype(int)  # shape (G, k)
        G, k = mat_sub.shape
        tot_subs = int(mat_sub.sum())
        mat_bin = (mat_sub > 0).astype(int)
        tot_bin = int(mat_bin.sum())

        if tot_subs == 0:
            continue

        def test_matrix(M: np.ndarray, total_val: int) -> tuple[float, float, float, float]:
            col_sums = M.sum(axis=0).astype(float)
            p = col_sums / total_val
            nz = p > 0
            obs_ent = -float(np.sum(p[nz] * np.log2(p[nz]))) / np.log2(k)
            obs_var = float(np.sum((col_sums - col_sums.mean()) ** 2))

            # Only permute rows that actually have non-constant entries (big speedup!)
            active_rows = [M[g_i] for g_i in range(G) if np.any(M[g_i] != M[g_i, 0])]
            const_sum = np.zeros(k, dtype=float)
            for g_i in range(G):
                if not np.any(M[g_i] != M[g_i, 0]):
                    const_sum += M[g_i]

            if not active_rows:
                return obs_ent, obs_var, 1.0, 1.0

            # Vectorized intra-row permutations in chunks
            n_act = len(active_rows)
            ge_var = 0
            le_ent = 0
            log2_k = np.log2(k)
            inv_tot = 1.0 / total_val

            # Generate random permutations per active row in batch (N_PERM, k)
            sim_col_sums = np.broadcast_to(const_sum, (N_PERM, k)).copy()
            for row_arr in active_rows:
                # Argsort random uniforms to get N_PERM independent permutations of 0..k-1
                idx = np.argsort(rng.random((N_PERM, k)), axis=1)
                sim_col_sums += row_arr[idx]

            sim_var = np.sum((sim_col_sums - (total_val / k)) ** 2, axis=1)
            ge_var = int(np.sum(sim_var >= obs_var - 1e-12))

            sim_p = sim_col_sums * inv_tot
            with np.errstate(divide="ignore", invalid="ignore"):
                log_p = np.where(sim_p > 0, np.log2(sim_p), 0.0)
            sim_ent = -np.sum(sim_p * log_p, axis=1) / log2_k
            le_ent = int(np.sum(sim_ent <= obs_ent + 1e-12))

            p_var = float((1 + ge_var) / (N_PERM + 1))
            p_ent = float((1 + le_ent) / (N_PERM + 1))
            return obs_ent, obs_var, p_var, p_ent

        ent_cnt, var_cnt, p_var_cnt, p_ent_cnt = test_matrix(mat_sub, tot_subs)
        ent_bin, var_bin, p_var_bin, p_ent_bin = test_matrix(mat_bin, tot_bin)

        rows.append({
            "KHIPU_ID": kid,
            "OKR_NUM": r.OKR_NUM,
            "kfg_name": r.kfg_name,
            "provenance": r.provenance,
            "det_family": r.det_family,
            "G": G,
            "k": k,
            "sub_in_table": tot_subs,
            "cords_with_sub": tot_bin,
            "col_sub_counts": str(mat_sub.sum(axis=0).tolist()),
            "col_bin_counts": str(mat_bin.sum(axis=0).tolist()),
            "ent_counts": round(ent_cnt, 4),
            "p_perm_counts_var": round(p_var_cnt, 5),
            "p_perm_counts_ent": round(p_ent_cnt, 5),
            "ent_binary": round(ent_bin, 4),
            "p_perm_binary_var": round(p_var_bin, 5),
            "p_perm_binary_ent": round(p_ent_bin, 5),
        })

    res = pd.DataFrame(rows)
    # Primary testable set: sub_in_table >= 5 (same n=38 as Step 2) + binary testable: cords_with_sub >= 5
    for mask_col, thresh, p_col, q_col in [
        ("sub_in_table", 5, "p_perm_counts_var", "q_perm_counts_var"),
        ("sub_in_table", 5, "p_perm_counts_ent", "q_perm_counts_ent"),
        ("cords_with_sub", 5, "p_perm_binary_var", "q_perm_binary_var"),
        ("cords_with_sub", 5, "p_perm_binary_ent", "q_perm_binary_ent"),
    ]:
        res[q_col] = np.nan
        m = res[mask_col] >= thresh
        if m.sum() > 0:
            res.loc[m, q_col] = false_discovery_control(res.loc[m, p_col].values, method="bh")

    res.to_csv(RESULTS_DIR / "step2_subsidiary_perm_5000.csv", index=False)

    print("\n=== PART B (#4): EXACT INTRA-ROW PERMUTATION TEST FOR SUBSIDIARY SLOTS (N_PERM=5,000) ===")
    for label, mask, p_col, q_col, ent_col in [
        ("Raw counts (sub_in_table>=5, var)",  res["sub_in_table"] >= 5,   "p_perm_counts_var", "q_perm_counts_var", "ent_counts"),
        ("Raw counts (sub_in_table>=5, ent)",  res["sub_in_table"] >= 5,   "p_perm_counts_ent", "q_perm_counts_ent", "ent_counts"),
        ("Binary presence (cords>=5, var)",    res["cords_with_sub"] >= 5, "p_perm_binary_var", "q_perm_binary_var", "ent_binary"),
        ("Binary presence (cords>=5, ent)",    res["cords_with_sub"] >= 5, "p_perm_binary_ent", "q_perm_binary_ent", "ent_binary"),
    ]:
        sub = res[mask]
        n_tot = len(sub)
        n_raw = int((sub[p_col] < 0.05).sum())
        n_fdr = int((sub[q_col] < 0.05).sum())
        # Corpus-level Fisher/Wilcoxon: excess concentration above uniform null (p < 0.5)
        n_lt_half = int((sub[p_col] < 0.5).sum())
        # Stouffer / Fisher combined p or binomial test on p < 0.05
        binom_p = stats.binomtest(n_raw, n_tot, 0.05, alternative="greater").pvalue if n_tot else np.nan
        print(
            f"  {label:38s}: n={n_tot:2d} | raw p<0.05: {n_raw:2d} ({100*n_raw/n_tot:4.1f}%) | "
            f"BH q<0.05: {n_fdr:2d} ({100*n_fdr/n_tot:4.1f}%) | med_ent={sub[ent_col].median():.3f} | "
            f"corpus binom p={binom_p:.2e}"
        )

    print("\n  Top khipus passing BH q<0.05 on BOTH raw counts AND binary presence:")
    both_sig = res[(res["q_perm_counts_var"] < 0.05) & (res["q_perm_binary_var"] < 0.05)].sort_values("p_perm_binary_var")
    if len(both_sig):
        print(
            both_sig[
                ["OKR_NUM", "kfg_name", "provenance", "G", "k", "sub_in_table", "cords_with_sub", "col_bin_counts", "q_perm_counts_var", "q_perm_binary_var"]
            ].round(4).to_string(index=False)
        )
    else:
        print("  (none passed both simultaneously)")
    return res


def compute_sz_calendar_candidates_okr(pendants: pd.DataFrame) -> set[int]:
    """
    Replicate the exact KFG S/Z spatial-block calendar criterion directly from OKR level-1 pendants:
      - pct_S >= 30% and pct_Z >= 30% (among all knots on level-1 pendants)
      - |s_moment - mid_moment| > 2 and |z_moment - mid_moment| > 2
      - strictly ordered: s_moment < mid_moment < z_moment OR z_moment < mid_moment < s_moment
    """
    cand_kids = set()
    for kid, g in pendants.groupby("KHIPU_ID", sort=False):
        g = g.reset_index(drop=True)
        pos = np.arange(1, len(g) + 1, dtype=float)
        ks = g["n_knots_S"].values.astype(float)
        kz = g["n_knots_Z"].values.astype(float)
        kt = g["n_knots"].values.astype(float)
        tot_k = kt.sum()
        if tot_k <= 0:
            continue
        pct_s = 100.0 * ks.sum() / tot_k
        pct_z = 100.0 * kz.sum() / tot_k
        if pct_s < 30.0 or pct_z < 30.0:
            continue
        s_mom = float(np.sum(pos * ks) / ks.sum())
        z_mom = float(np.sum(pos * kz) / kz.sum())
        mid_mom = float(np.sum(pos * kt) / tot_k)
        if abs(s_mom - mid_mom) > 2.0 and abs(z_mom - mid_mom) > 2.0:
            if (s_mom < mid_mom < z_mom) or (z_mom < mid_mom < s_mom):
                cand_kids.add(int(kid))
    return cand_kids


def run_part_c_outlier_enrichment(pendants: pd.DataFrame, sig86: pd.DataFrame) -> pd.DataFrame:
    """Exact Fisher enrichment tests: Deterministic OUT (21/86) and HDBSCAN noise (45/86) vs special khipus."""
    sz_cal_kids = compute_sz_calendar_candidates_okr(pendants)

    df = sig86.copy()
    df["is_official_cal10"] = df["OKR_NUM"].isin(OFFICIAL_CALENDAR_KFG10)
    df["is_official_cal11"] = df["OKR_NUM"].isin(OFFICIAL_CALENDAR_PLUS_UR006)
    df["is_sz_cal_cand"] = df["KHIPU_ID"].isin(sz_cal_kids)
    df["is_partb_new_pair"] = df["kfg_name"].isin(PARTB_NEW_PAIR_NAMES) | df["OKR_NUM"].isin(PARTB_NEW_PAIR_NAMES)
    df["is_any_special"] = df["is_official_cal11"] | df["is_sz_cal_cand"] | df["is_partb_new_pair"]

    rows = []
    for grp_name, grp_mask in [
        ("Deterministic OUT (21/86)", df["is_det_outlier"].values),
        ("HDBSCAN noise -1 (45/86)", df["is_hdbscan_noise"].values),
        ("Both OUT & HDBSCAN -1 (9/86)", (df["is_det_outlier"] & df["is_hdbscan_noise"]).values),
    ]:
        for target_name, t_col in [
            ("Official KFG calendar (10 list)", "is_official_cal10"),
            ("Official calendar + UR006/KH0242", "is_official_cal11"),
            ("S/Z spatial-block calendar cand", "is_sz_cal_cand"),
            ("KFG Part B new paired khipus (9 pairs)", "is_partb_new_pair"),
            ("UNION (any calendar OR new pair)", "is_any_special"),
        ]:
            t = df[t_col].values
            a = int(np.sum(grp_mask & t))
            b = int(np.sum(grp_mask & (~t)))
            c = int(np.sum((~grp_mask) & t))
            d = int(np.sum((~grp_mask) & (~t)))
            or_val, p_two = stats.fisher_exact([[a, b], [c, d]], alternative="two-sided")
            _, p_gt = stats.fisher_exact([[a, b], [c, d]], alternative="greater")
            rows.append({
                "outlier_def": grp_name,
                "target": target_name,
                "in_outliers": f"{a}/{a+b}",
                "in_inliers": f"{c}/{c+d}",
                "pct_outliers": round(100.0 * a / max(a + b, 1), 1),
                "pct_inliers": round(100.0 * c / max(c + d, 1), 1),
                "odds_ratio": round(float(or_val), 2) if np.isfinite(or_val) else np.inf,
                "fisher_p_two": round(float(p_two), 5),
                "fisher_p_greater": round(float(p_gt), 5),
            })

    out_df = pd.DataFrame(rows)
    out_df.to_csv(RESULTS_DIR / "step2_outlier_enrichment_fisher.csv", index=False)

    print(f"\n=== PART C (#5): FISHER EXACT ENRICHMENT OF OUTLIERS (OKR S/Z cal total={len(sz_cal_kids)}) ===")
    print(out_df.to_string(index=False))

    print("\n  Exact list of special khipus inside the 86 doubly-tabular cohort:")
    sp = df[df["is_any_special"]][
        ["OKR_NUM", "kfg_name", "provenance", "G", "k", "det_family", "hdbscan_cluster",
         "is_official_cal11", "is_sz_cal_cand", "is_partb_new_pair"]
    ]
    print(sp.to_string(index=False))
    return out_df


def main():
    pendants, _ = load_okr_pendants()
    cohorts = load_cohorts()
    sig86 = pd.read_csv(RESULTS_DIR / "step2_signatures_86.csv")

    run_part_a_step1_audit(pendants, cohorts)
    run_part_b_subsidiary_perm(pendants, sig86)
    run_part_c_outlier_enrichment(pendants, sig86)


if __name__ == "__main__":
    main()
