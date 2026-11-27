"""
gate0_santa_valley.py — Pre-registered Gate 0 (Disjoint Positive Control):
  - Gate 0A: Exact computational replication of Medrano & Urton (2018, Ethnohistory 65:1)
    directly from OKR v2.1.0 + KFG Excel workbooks (UR087-UR092 / KH0323-KH0328):
      * 133 six-cord color-banded groups (Table 1: 63 Recto / 70 Verso first cords)
      * Sum of first-cord numeric values (Note 17: 359 observed -> 367/368 with broken cords)
      * Average first-cord tribute values by khipu & moiety (Figures 3a & 3b)
      * First-cord distinct color count & singletons (pp. 18-19: 32 colors, 17 singletons
        vs 30 distinct first names, 14 singletons in Table 4).
  - Gate 0B: Calibration of the frozen Step 3 alignment metric A_JSD + Null 1 (exact k=6,
    720 permutations) + Null 2 (corpus-empirical) on within-Santa-Valley cross-khipu match
    (UR087/KH0323 34x6 vs UR089/KH0325 4x6, strictly disjoint from the 86 doubly-tabular set).
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from src.data_loader import RESULTS_DIR, extract_modal_table, load_cohorts, load_kfg_excel_pendants, load_okr_pendants
from src.step3_colonial_alignment import align_jsd, align_spearman, mean_column_shares, null1_column_perm_p, null2_empirical_p

SANTA_MAP = [
    ("UR087", "KH0323", 48, 0),
    ("UR092", "KH0328", 10, 0),
    ("UR088", "KH0324", 0, 9),
    ("UR089", "KH0325", 0, 34),
    ("UR091", "KH0327", 0, 15),
    ("UR090", "KH0326", 5, 12),
]


def segment_six_cord_groups(df_okr: pd.DataFrame, df_kfg: pd.DataFrame, kfg_name: str) -> pd.DataFrame:
    """
    Partition a Santa Valley khipu's pendants into its Medrano & Urton (2018) six-cord groups.
    Per Note 9 (p. 21): 804 pendants total; on UR090 cords 31 & 38 are excluded and cords 57-66
    form one 10-cord group (so UR090 has 17 groups instead of 18).
    Aligns OKR attachment/color with KFG vertical-position values by pendant index (1..N).
    """
    n = min(len(df_okr), len(df_kfg))
    merged = pd.DataFrame({
        "p_idx": np.arange(1, n + 1),
        "okr_val": df_okr["cord_val"].values[:n].astype(int),
        "kfg_val": df_kfg["kfg_val"].values[:n].astype(int),
        "att": df_okr["att_type"].values[:n],
        "color": df_okr["color_norm"].values[:n],
        "kfg_term": df_kfg["kfg_term"].values[:n],
    })

    # Identify first cord of each six-cord group (1-indexed)
    if kfg_name == "UR090":
        # Per Note 9 (p. 21): cords 31 and 38 excluded; cords 57-66 are a single 10-cord group
        # Let's segment via color runs / exact Note 9 offsets:
        # Cords 1..30 -> 5 groups (1,7,13,19,25); skip 31; 32..37 -> 1 group (32); skip 38;
        # 39..56 -> 3 groups (39,45,51); 57..66 -> 1 ten-cord group (57);
        # 67..108 -> 7 groups (67,73,79,85,91,97,103). Total = 5+1+3+1+7 = 17 groups!
        first_idxs = [1, 7, 13, 19, 25, 32, 39, 45, 51, 57, 67, 73, 79, 85, 91, 97, 103]
    elif kfg_name == "UR087":
        # 48 six-cord groups (= 288 cords; remaining pendants are non-group/preamble/postamble)
        # Use exact contiguous 6-cord blocks starting from first group in KFG clusters of size 6 (or 1..288)
        cl_sz = df_kfg.groupby("kfg_cluster").size()
        # Check where first size-6 cluster begins
        six_cids = cl_sz[cl_sz >= 6].index.tolist()
        first_p = int(df_kfg[df_kfg["kfg_cluster"] == six_cids[0]]["p_idx"].min()) if six_cids else 1
        # Collect all 6-cord starts across clusters whose size is a multiple of 6
        first_idxs = []
        for cid, grp in df_kfg[df_kfg["kfg_cluster"] >= 0].groupby("kfg_cluster", sort=False):
            sz = len(grp)
            if sz >= 6 and sz % 6 == 0:
                p0 = int(grp["p_idx"].iloc[0])
                for off in range(0, sz, 6):
                    first_idxs.append(p0 + off)
    elif kfg_name == "UR089":
        # 34 six-cord groups (= 204 cords; 206 total pendants, 2 singletons)
        first_idxs = []
        for cid, grp in df_kfg[df_kfg["kfg_cluster"] >= 0].groupby("kfg_cluster", sort=False):
            sz = len(grp)
            if sz >= 6 and sz % 6 == 0:
                p0 = int(grp["p_idx"].iloc[0])
                for off in range(0, sz, 6):
                    first_idxs.append(p0 + off)
    else:
        # UR088 (54=9x6), UR091 (90=15x6), UR092 (60=10x6)
        first_idxs = list(range(1, n + 1, 6))

    sub = merged[merged["p_idx"].isin(first_idxs)].copy()
    sub["kfg_name"] = kfg_name
    return sub


def main():
    pendants, _ = load_okr_pendants()
    cohorts = load_cohorts()
    pend_by_okr = {str(okr): g.reset_index(drop=True) for okr, g in pendants.groupby("OKR_NUM", sort=False)}

    print("=== GATE 0A: REPLICATION OF MEDRANO & URTON (2018) TABLE 1, NOTE 17 & COLOR COUNTS ===")
    first_cords_all = []
    summary_rows = []
    for kfg_name, okr_num, pub_R, pub_V in SANTA_MAP:
        df_okr = pend_by_okr[okr_num]
        df_kfg = load_kfg_excel_pendants(kfg_name)
        fc = segment_six_cord_groups(df_okr, df_kfg, kfg_name)
        first_cords_all.append(fc)

        n_R = int((fc["att"] == "R").sum())
        n_V = int((fc["att"] == "V").sum())
        n_U = int((~fc["att"].isin(["R", "V"])).sum())
        sum_kfg = int(fc["kfg_val"].sum())
        sum_okr = int(fc["okr_val"].sum())
        n_broken = int((fc["kfg_term"] == "B").sum())
        mean_kfg = float(fc["kfg_val"].mean())
        summary_rows.append({
            "kfg_name": kfg_name,
            "OKR_NUM": okr_num,
            "n_groups_obs": len(fc),
            "n_groups_pub": pub_R + pub_V,
            "R_obs": n_R,
            "R_pub": pub_R,
            "V_obs": n_V,
            "V_pub": pub_V,
            "U_obs": n_U,
            "first_cord_sum_kfg": sum_kfg,
            "first_cord_sum_okr": sum_okr,
            "first_cord_mean_kfg": round(mean_kfg, 2),
            "n_broken_first_cords": n_broken,
        })

    sum_df = pd.DataFrame(summary_rows)
    all_fc = pd.concat(first_cords_all, ignore_index=True)
    print(sum_df.to_string(index=False))

    tot_groups = int(sum_df["n_groups_obs"].sum())
    tot_R = int(sum_df["R_obs"].sum())
    tot_V = int(sum_df["V_obs"].sum())
    tot_sum_kfg = int(sum_df["first_cord_sum_kfg"].sum())
    tot_sum_okr = int(sum_df["first_cord_sum_okr"].sum())
    col_vc = all_fc["color"].value_counts()
    n_uniq_col = len(col_vc)
    n_single_col = int((col_vc == 1).sum())

    print(f"\nTotals across all 6 Santa Valley khipus:")
    print(f"  Six-cord groups:     {tot_groups} observed vs 133 published")
    print(f"  First-cord R / V:    {tot_R} R / {tot_V} V observed vs 63 R / 70 V published")
    print(f"  First-cord val sum:  {tot_sum_kfg} (KFG Excel) / {tot_sum_okr} (OKR base-10) vs 359 observed (367/368 w/ broken) in Note 17")
    print(f"  Distinct 1st colors: {n_uniq_col} ({n_single_col} singletons) vs 32 (17 singletons) published (pp. 18-19)")
    sum_df.to_csv(RESULTS_DIR / "gate0a_santa_valley_replication.csv", index=False)

    # --- GATE 0B: Calibration of frozen Step 3 A_JSD + Null 1 (exact k=6) + Null 2 ---
    print("\n=== GATE 0B: FROZEN STEP 3 ALIGNMENT CALIBRATION ON SANTA VALLEY (k=6, DISJOINT FROM 86) ===")
    # Extract exact 6-column matrices from KFG Excel for UR087 (34x6) and UR089 (4x6), plus ALL 133 six-cord groups pooled!
    def kfg_k6_matrix(kfg_name):
        df = load_kfg_excel_pendants(kfg_name)
        mats = []
        for _, grp in df[df["kfg_cluster"] >= 0].groupby("kfg_cluster", sort=False):
            sz = len(grp)
            if sz >= 6 and sz % 6 == 0:
                vals = grp["kfg_val"].values.astype(float)
                for off in range(0, sz, 6):
                    mats.append(vals[off:off+6])
        return np.array(mats, dtype=float)

    m87 = kfg_k6_matrix("UR087")
    m89 = kfg_k6_matrix("UR089")
    m91 = kfg_k6_matrix("UR091")
    m92 = kfg_k6_matrix("UR092")
    m88 = kfg_k6_matrix("UR088")

    # Build empirical pool of all k=6 khipus in OKR (outside Santa Valley) for Null 2
    emp_pool_k6 = []
    for kid, g in pendants.groupby("KHIPU_ID", sort=False):
        if g["OKR_NUM"].iloc[0] in {"KH0323", "KH0324", "KH0325", "KH0326", "KH0327", "KH0328"}:
            continue
        if len(g) < 10:
            continue
        tab = extract_modal_table(g)
        if tab is not None and tab["k"] == 6:
            emp_pool_k6.append(tab["mat_val"])
    print(f"Empirical non-Santa k=6 khipu pool: {len(emp_pool_k6)} khipus")

    rng = np.random.default_rng(1670)
    cal_rows = []
    sv_mats = {"UR087": m87, "UR088": m88, "UR089": m89, "UR091": m91, "UR092": m92}
    for name_a, name_b in [("UR087", "UR089"), ("UR087", "UR091"), ("UR089", "UR091"), ("UR087", "UR088"), ("UR087", "UR092")]:
        ma, mb = sv_mats[name_a], sv_mats[name_b]
        pa, pb = mean_column_shares(ma), mean_column_shares(mb)
        a_jsd, orient = align_jsd(pa, pb)
        a_rho = align_spearman(pa, pb)
        p1 = null1_column_perm_p(pa, pb, a_jsd, rng)
        p2 = null2_empirical_p(ma, pb, a_jsd, emp_pool_k6, rng)
        cal_rows.append({
            "khipu_a": name_a,
            "khipu_b": name_b,
            "rows_a": len(ma),
            "rows_b": len(mb),
            "shares_a": np.round(pa, 3).tolist(),
            "shares_b": np.round(pb, 3).tolist(),
            "A_JSD": round(a_jsd, 4),
            "orient": orient,
            "A_rank": round(a_rho, 4),
            "p_null1_exact": round(p1, 5),
            "p_null2_emp": round(p2, 5),
            "p_both": round(max(p1, p2), 5),
        })
    cal_df = pd.DataFrame(cal_rows)
    cal_df.to_csv(RESULTS_DIR / "gate0b_santa_valley_calibration.csv", index=False)
    print(cal_df[["khipu_a", "khipu_b", "rows_a", "rows_b", "A_JSD", "orient", "A_rank", "p_null1_exact", "p_null2_emp", "p_both"]].to_string(index=False))


if __name__ == "__main__":
    main()
