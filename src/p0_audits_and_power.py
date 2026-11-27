"""
p0_audits_and_power.py — P0 Credibility & Robustness Audits:
  1. Pre-registered secondary Step 3 analysis for |k_K - k_D| == 1 (with unbiased max-slice
     Null 1 and Null 2, N_PERM = 5,000, BH-FDR q < 0.05).
  2. Statistical power curve for Step 3 (A_JSD recovery under increasing relative noise sigma
     on real k=5, 6, 7 khipus under Null 1, Null 2, and BH-FDR thresholds).
  3. Deduplicated Phase 2A re-run (excluding Berlin, Huacho, Nazca, Huaquerones, Puruchuco,
     and Incahuasi duplicate/paired specimens) with standardized effect sizes, 95% CIs,
     exact permutation bounds (p <= 1/(N_PERM+1)), and Holm/BH correction across all 5
     candidate subsidiary regularities.
  4. Explicit Part 2 (Groups 22-24, cords 105-119) verification of UR1175 (VA 42518) against
     all three khipus in Berlin lot VA 63038 (UR231/KH0468, UR232/KH0469, UR233/KH0470).
"""

from __future__ import annotations

import itertools
import json
import sqlite3
import warnings

import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import binomtest, false_discovery_control, norm, spearmanr, wilcoxon

from src.data_loader import (
    DATA_DIR,
    RESULTS_DIR,
    extract_modal_table,
    find_okr_db,
    load_cohorts,
    load_okr_pendants,
)
from src.phase2a_subsidiary_grammar import (
    build_okr_modal_tables,
    load_okr_cords_and_subs,
    run_law1_slot_presence,
    run_law2_4th_axis_color,
    run_law3_5th_axis_value,
    run_law5_doseresponse,
    wilson_ci,
)
from src.phase2c_pair_discovery import load_kfg_pendants
from src.step3_colonial_alignment import (
    align_jsd_batch,
    align_spearman,
    mean_column_shares,
)

warnings.filterwarnings("ignore")

SCHEMA_DIR = DATA_DIR / "colonial_schemas"
SEED = 20261127

# Explicit deduplication drop list (keeping one representative per known/putative duplicate or lot pair):
#   - KH0483 (UR247, dup of KH0140/AS125)
#   - KH0304 (UR067, dup of KH0303/UR066)
#   - KH0473 (UR237, reversed dup of KH0121/AS107)
#   - KH0469 (UR232) & KH0470 (UR233) (putative Huacho/VA 63038 companions of KH0192/UR1175)
#   - KH0156 (UR1140, consecutive Berlin Nazca museum number VA16636 paired with KH0152/VA16635)
#   - KH0305 (UR068, Puruchuco companion of KH0301/UR064)
#   - KH0503 (UR267A) & KH0504 (UR267B) (Incahuasi pairs of KH0491/UR255 & KH0492/UR256)
DEDUP_DROP_OKR = {
    "KH0483", "KH0304", "KH0473", "KH0469", "KH0470",
    "KH0156", "KH0305", "KH0503", "KH0504",
}


def bootstrap_median_ci(arr: np.ndarray, rng: np.random.Generator, n_boot: int = 2000) -> tuple[float, float]:
    arr = np.asarray(arr, dtype=float)
    arr = arr[~np.isnan(arr)]
    if len(arr) == 0:
        return (np.nan, np.nan)
    idx = rng.integers(0, len(arr), size=(n_boot, len(arr)))
    meds = np.median(arr[idx], axis=1)
    return (float(np.quantile(meds, 0.025)), float(np.quantile(meds, 0.975)))


def wilcoxon_effect_and_z(diff: np.ndarray) -> tuple[float, float, float]:
    """Return (z_stat, rank_biserial_r, p_asymp_or_exact) cleanly without p=0 underflow."""
    d = np.asarray(diff, dtype=float)
    d = d[~np.isnan(d)]
    d = d[d != 0]
    n = len(d)
    if n < 5:
        return (np.nan, np.nan, np.nan)
    res = wilcoxon(d, method="approx" if n >= 50 else "exact")
    # Compute rank-biserial correlation r_rb = (W+ - W-) / (n*(n+1)/2)
    ranks = stats.rankdata(np.abs(d))
    w_pos = float(ranks[d > 0].sum())
    w_neg = float(ranks[d < 0].sum())
    tot = n * (n + 1) / 2.0
    r_rb = (w_pos - w_neg) / tot
    # Standardized Z (with tie/variance formula)
    sigma = np.sqrt(n * (n + 1) * (2 * n + 1) / 24.0)
    z = (w_pos - tot / 2.0) / sigma if sigma > 0 else np.nan
    p = float(res.pvalue)
    return (float(z), float(r_rb), p)


# =====================================================================
# 1. STEP 3 SECONDARY SENSITIVITY: |k_K - k_D| == 1
# =====================================================================
def _align_pm1_batch(P_k: np.ndarray, P_d: np.ndarray) -> np.ndarray:
    """
    Unbiased |k_K - k_D| == 1 alignment:
    Drop either the first or last column of the wider table, re-normalize the resulting
    k_min-column share vector, and take max A_JSD across the two boundary drops.
    Works for P_k (N, k_k) and P_d (N, k_d) or (1, k_d).
    """
    P_k = np.asarray(P_k, dtype=float)
    P_d = np.asarray(P_d, dtype=float)
    if P_k.ndim == 1:
        P_k = P_k[None, :]
    if P_d.ndim == 1:
        P_d = P_d[None, :]
    kk, kd = P_k.shape[1], P_d.shape[1]
    if kk == kd + 1:
        s1 = align_jsd_batch(P_k[:, 1:], P_d)
        s2 = align_jsd_batch(P_k[:, :-1], P_d)
        return np.maximum(s1, s2)
    elif kd == kk + 1:
        s1 = align_jsd_batch(P_k, P_d[:, 1:])
        s2 = align_jsd_batch(P_k, P_d[:, :-1])
        return np.maximum(s1, s2)
    raise ValueError(f"Expected |k_k - k_d| == 1, got {kk} vs {kd}")


def run_step3_secondary_pm1(
    cohorts: pd.DataFrame, okr_pends: pd.DataFrame, rng: np.random.Generator, n_perm: int = 5000
) -> pd.DataFrame:
    schemas = []
    for p in sorted(SCHEMA_DIR.glob("*.json")):
        with open(p, encoding="utf-8") as f:
            schemas.append(json.load(f))

    # Build empirical row pools per k for Null 2
    emp_pool: dict[int, list[np.ndarray]] = {}
    khipu_tables: dict[int, np.ndarray] = {}
    by_kid = {int(kid): g for kid, g in okr_pends.groupby("KHIPU_ID", sort=False)}
    for row in cohorts.itertuples(index=False):
        if pd.isna(row.h3_m) or int(row.h3_m) < 3:
            continue
        kid = int(row.KHIPU_ID)
        if kid not in by_kid:
            continue
        tbl = extract_modal_table(by_kid[kid])
        if tbl is None or tbl["G"] < 2:
            continue
        k = int(tbl["k"])
        mat = tbl["mat_val"]
        khipu_tables[kid] = mat
        rs = mat.sum(axis=1)
        nz = rs > 0
        if np.any(nz):
            normed = mat[nz] / rs[nz, None]
            emp_pool.setdefault(k, []).extend(normed)

    rows = []
    for cohort_name, mask_col in [("primary_86", "doubly_tabular"), ("sens_H3b_only", "h3b_only")]:
        sub_coh = cohorts[cohorts[mask_col]].copy()
        for sc in schemas:
            mat_d = np.asarray(sc["value_matrix"], dtype=float)
            kd = int(mat_d.shape[1])
            p_d = mean_column_shares(mat_d)
            for row in sub_coh.itertuples(index=False):
                if pd.isna(row.h3_m):
                    continue
                kk = int(row.h3_m)
                if abs(kk - kd) != 1:
                    continue
                kid = int(row.KHIPU_ID)
                if kid not in khipu_tables:
                    continue
                mat_k = khipu_tables[kid]
                p_k = mean_column_shares(mat_k)
                if np.any(np.isnan(p_k)):
                    continue

                obs = float(_align_pm1_batch(p_k, p_d)[0])

                # Unbiased Null 1: permute columns of D (using same max-boundary-drop operator!)
                if kd <= 7:
                    perms = np.array(list(itertools.permutations(range(kd))), dtype=int)
                    scores1 = _align_pm1_batch(p_k[None, :], p_d[perms])
                    p_n1 = float(np.mean(scores1 >= obs - 1e-12))
                else:
                    idx = np.argsort(rng.random((n_perm, kd)), axis=1)
                    scores1 = _align_pm1_batch(p_k[None, :], p_d[idx])
                    p_n1 = float((1 + int(np.sum(scores1 >= obs - 1e-12))) / (n_perm + 1))

                # Unbiased Null 2: corpus-empirical pool if >=20 khipus share kk, else within-khipu row shuffles
                pool = emp_pool.get(kk, [])
                rs = mat_k.sum(axis=1)
                nz = rs > 0
                norm_k = mat_k[nz] / rs[nz, None]
                G_nz = norm_k.shape[0]
                if len(pool) >= 20 * 3:
                    pool_arr = np.asarray(pool, dtype=float)
                    draw_idx = rng.integers(0, len(pool_arr), size=(n_perm, G_nz))
                    syn = pool_arr[draw_idx].mean(axis=1)
                else:
                    # within-khipu row-wise column shuffles
                    syn = np.empty((n_perm, kk), dtype=float)
                    for b in range(n_perm):
                        perm_idx = np.argsort(rng.random((G_nz, kk)), axis=1)
                        shuf = np.take_along_axis(norm_k, perm_idx, axis=1)
                        syn[b] = shuf.mean(axis=0)
                scores2 = _align_pm1_batch(syn, p_d[None, :])
                p_n2 = float((1 + int(np.sum(scores2 >= obs - 1e-12))) / (n_perm + 1))

                k_min = min(kk, kd)
                p_both = max(p_n1, p_n2) if k_min >= 6 else p_n2
                rows.append({
                    "cohort": cohort_name,
                    "OKR_NUM": row.OKR_NUM,
                    "provenance": row.provenance_clean,
                    "k_khipu": kk,
                    "G_khipu": mat_k.shape[0],
                    "schema_id": sc["schema_id"],
                    "k_schema": kd,
                    "k_min": k_min,
                    "A_JSD_pm1": round(obs, 5),
                    "p_null1": p_n1,
                    "p_null2": p_n2,
                    "p_eval": p_both,
                })

    df = pd.DataFrame(rows)
    df["q_bh"] = np.nan
    for coh in df["cohort"].unique():
        m = df["cohort"] == coh
        df.loc[m, "q_bh"] = false_discovery_control(df.loc[m, "p_eval"].values, method="bh")
    df["pass_fdr"] = df["q_bh"] < 0.05
    return df


# =====================================================================
# 2. STEP 3 POWER ANALYSIS (NOISE PERTURBATION CURVE FOR 0/30)
# =====================================================================
def run_step3_power_analysis(
    cohorts: pd.DataFrame, okr_pends: pd.DataFrame, rng: np.random.Generator, n_trials: int = 120, n_perm_n2: int = 1000
) -> pd.DataFrame:
    """
    Perturb real tabular khipu column share vectors (for k=5, 6, 7) with multiplicative log-normal
    noise of increasing relative SD sigma in [0.00..1.00], and measure empirical detection power
    under Null 1, Null 2, joint intersection-union (p_both < 0.05), and Bonferroni/FDR proxy
    (p < 0.05/30 = 0.00167 for m=30 tests).
    """
    by_kid = {int(kid): g for kid, g in okr_pends.groupby("KHIPU_ID", sort=False)}
    rep_tables: dict[int, list[tuple[str, np.ndarray]]] = {5: [], 6: [], 7: []}
    emp_pool: dict[int, list[np.ndarray]] = {5: [], 6: [], 7: []}
    for row in cohorts.itertuples(index=False):
        if pd.isna(row.h3_m):
            continue
        kid = int(row.KHIPU_ID)
        if kid not in by_kid:
            continue
        tbl = extract_modal_table(by_kid[kid])
        if tbl is None or tbl["G"] < 4 or int(tbl["k"]) not in (5, 6, 7):
            continue
        k = int(tbl["k"])
        mat = tbl["mat_val"]
        rs = mat.sum(axis=1)
        nz = rs > 0
        if nz.sum() < 4:
            continue
        normed = mat[nz] / rs[nz, None]
        emp_pool[k].extend(normed)
        if row.doubly_tabular and len(rep_tables[k]) < 4:
            p_k = normed.mean(axis=0)
            if np.std(p_k) > 0.04:  # non-flat real tabular profile
                rep_tables[k].append((row.OKR_NUM, normed))

    sigmas = [0.0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.40, 0.50, 0.75, 1.00]
    rows = []
    for k in (5, 6, 7):
        perms = np.array(list(itertools.permutations(range(k))), dtype=int)
        pool_arr = np.asarray(emp_pool[k], dtype=float)
        khipu_list = rep_tables[k]
        for sigma in sigmas:
            a_jsds, p1s, p2s, pbs = [], [], [], []
            for trial in range(n_trials):
                okr_num, norm_k = khipu_list[trial % len(khipu_list)]
                p_k = norm_k.mean(axis=0)
                # Multiplicative log-normal noise on column shares (re-normalized to simplex)
                if sigma == 0.0:
                    p_d = p_k.copy()
                else:
                    mult = np.exp(rng.normal(0.0, sigma, size=k))
                    p_d = (p_k * mult)
                    p_d = p_d / p_d.sum()

                obs = float(align_jsd_batch(p_k[None, :], p_d[None, :])[0])
                # Exact Null 1
                s1 = align_jsd_batch(p_k[None, :], p_d[perms])
                p1 = float(np.mean(s1 >= obs - 1e-12))
                # Empirical Null 2
                G_nz = norm_k.shape[0]
                draw_idx = rng.integers(0, len(pool_arr), size=(n_perm_n2, G_nz))
                syn = pool_arr[draw_idx].mean(axis=1)
                s2 = align_jsd_batch(syn, p_d[None, :])
                p2 = float((1 + int(np.sum(s2 >= obs - 1e-12))) / (n_perm_n2 + 1))

                a_jsds.append(obs)
                p1s.append(p1)
                p2s.append(p2)
                pbs.append(max(p1, p2) if k >= 6 else p2)

            a_jsds = np.array(a_jsds); p1s = np.array(p1s); p2s = np.array(p2s); pbs = np.array(pbs)
            rows.append({
                "k": k,
                "n_rep_khipus": len(khipu_list),
                "rep_khipus": ",".join(x[0] for x in khipu_list),
                "noise_sigma": sigma,
                "n_trials": n_trials,
                "mean_A_JSD": round(float(a_jsds.mean()), 4),
                "sd_A_JSD": round(float(a_jsds.std()), 4),
                "power_null1_p05": round(float((p1s < 0.05).mean()), 3),
                "power_null2_p05": round(float((p2s < 0.05).mean()), 3),
                "power_eval_p05": round(float((pbs < 0.05).mean()), 3),
                "power_eval_p005": round(float((pbs < 0.005).mean()), 3),
                "power_fdr30_proxy": round(float((pbs < (0.05 / 30.0)).mean()), 3),
                "min_possible_null1_p": round(2.0 / len(perms), 5),
            })
    return pd.DataFrame(rows)


# =====================================================================
# 3. DEDUPLICATED PHASE 2A RE-RUN + EFFECT SIZES + MULTI-TEST FDR
# =====================================================================
def run_phase2a_dedup_audit(rng: np.random.Generator) -> tuple[pd.DataFrame, pd.DataFrame]:
    okr_m, cords, subs = load_okr_cords_and_subs()
    cords_ded = cords[~cords["OKR_NUM"].isin(DEDUP_DROP_OKR)].copy()
    subs_ded = subs[~subs["OKR_NUM"].isin(DEDUP_DROP_OKR)].copy()
    p1_mod_ded, sub_mod_ded = build_okr_modal_tables(cords_ded, subs_ded)

    # Empty KFG dict to run OKR-dedup cleanly
    empty_kdf = pd.DataFrame(columns=["kfg_name", "level", "parent_name", "color", "val", "cord_name"])

    # R1: Slot-locked presence
    l1 = run_law1_slot_presence(p1_mod_ded, empty_kdf, {}, rng, n_perm=2000)
    n1_sig = int((l1["q_bh"] < 0.05).sum())
    n1_tot = len(l1)
    p1_hat, p1_lo, p1_hi = wilson_ci(n1_sig, n1_tot)
    p1_bin = float(binomtest(n1_sig, n1_tot, 0.05, alternative="greater").pvalue)

    # R2: 4th-axis subsidiary color (strict sub != parent)
    l2 = run_law2_4th_axis_color(sub_mod_ded, empty_kdf, {}, rng)
    ind2 = l2["delta_nonparent"].dropna().values
    z2, r_rb2, p2_w = wilcoxon_effect_and_z(ind2)
    ci2_lo, ci2_hi = bootstrap_median_ci(ind2 * 100.0, rng)

    # R3: 5th-axis subsidiary value homology
    l3 = run_law3_5th_axis_value(sub_mod_ded, empty_kdf, {}, rng, n_perm=1000)
    exc3 = l3["excess_r2"].values
    z3, r_rb3, p3_w = wilcoxon_effect_and_z(exc3)
    ci3_lo, ci3_hi = bootstrap_median_ci(exc3 * 100.0, rng)

    # R4a: L1->L2 hierarchical decay (pair-level AND khipu-level)
    g4 = subs_ded[(subs_ded["CORD_LEVEL"] == 2) & (subs_ded["p_CORD_LEVEL"] == 1) & (subs_ded["cord_val"] > 0) & (subs_ded["p_cord_val"] > 0)]
    z4p, r_rb4p, p4p = wilcoxon_effect_and_z(g4["cord_val"].values - g4["p_cord_val"].values)
    kh4 = g4.groupby("KHIPU_ID").apply(lambda x: (x["cord_val"] < x["p_cord_val"]).mean()).values
    z4k, r_rb4k, p4k = wilcoxon_effect_and_z(kh4 - 0.5)
    ci4k_lo, ci4k_hi = bootstrap_median_ci(kh4 * 100.0, rng)

    # R4c: Non-trivial exact part-whole (parent >= 10, >=2 knotted subs), N_PERM = 5,000
    kn_ded = subs_ded[subs_ded["cord_val"] > 0]
    ag = kn_ded.groupby(["KHIPU_ID", "OKR_NUM", "direct_parent", "p_cord_val"]).agg(
        n_ksubs=("CORD_ID", "count"), sum_ksubs=("cord_val", "sum")
    ).reset_index()
    m10 = ag[(ag["n_ksubs"] >= 2) & (ag["p_cord_val"] >= 10)]
    obs_pw = int((m10["sum_ksubs"] == m10["p_cord_val"]).sum())
    pw_groups = [(grp["sum_ksubs"].values, grp["p_cord_val"].values) for _, grp in m10.groupby("KHIPU_ID") if len(grp) >= 2]
    n_perm_pw = 5000
    null_pw = np.empty(n_perm_pw, dtype=int)
    for it in range(n_perm_pw):
        tot = 0
        for s_arr, p_arr in pw_groups:
            tot += int((s_arr == rng.permutation(p_arr)).sum())
        null_pw[it] = tot
    p_pw = float((int(np.sum(null_pw >= obs_pw)) + 1) / (n_perm_pw + 1))

    # R5: Vertical descending sibling regularity (dose-response dAP >= 3cm vs 0cm neg control + khipu-level)
    dr_ded, cat_ded, meta5_ded = run_law5_doseresponse(subs_ded)
    r5_0 = dr_ded.iloc[0]
    r5_3 = dr_ded.iloc[3]

    # Build deduplicated effect-size table & apply BH-FDR + Holm across the 5 primary regularity tests
    reg_rows = [
        {
            "regularity": "R1: Slot-locked subsidiary presence (BH q<0.05)",
            "sample_dedup": f"{n1_sig}/{n1_tot} khipus",
            "effect_metric": "Proportion BH-significant",
            "effect_size": f"{100*p1_hat:.1f}% (95% CI {100*p1_lo:.1f}-{100*p1_hi:.1f}%)",
            "test_type": "Exact binomial vs 5% nominal",
            "raw_p": p1_bin,
            "p_display": f"{p1_bin:.2e}",
        },
        {
            "regularity": "R2: 4th-axis subsidiary color (strict sub!=parent)",
            "sample_dedup": f"{int((ind2>0).sum())}/{len(ind2)} khipus",
            "effect_metric": "Median delta (pp); rank-biserial r",
            "effect_size": f"+{np.median(ind2)*100:.2f} pp (95% CI {ci2_lo:+.2f} to {ci2_hi:+.2f}); r_rb={r_rb2:+.3f} (Z={z2:+.2f})",
            "test_type": "Wilcoxon signed-rank (asymptotic)",
            "raw_p": p2_w,
            "p_display": f"{p2_w:.2e}",
        },
        {
            "regularity": "R3: 5th-axis subsidiary numeric homology",
            "sample_dedup": f"{int((exc3>0).sum())}/{len(exc3)} khipus",
            "effect_metric": "Median excess R^2 (pp); rank-biserial r",
            "effect_size": f"+{np.median(exc3)*100:.2f} pp (95% CI {ci3_lo:+.2f} to {ci3_hi:+.2f}); r_rb={r_rb3:+.3f} (Z={z3:+.2f})",
            "test_type": "Wilcoxon signed-rank (asymptotic)",
            "raw_p": p3_w,
            "p_display": f"{p3_w:.2e}",
        },
        {
            "regularity": "R4: Hierarchical decay L1->L2 (khipu-level)",
            "sample_dedup": f"{len(kh4)} khipus ({len(g4)} pairs)",
            "effect_metric": "Median khipu frac(L2<L1); rank-biserial r",
            "effect_size": f"{np.median(kh4)*100:.1f}% (95% CI {ci4k_lo:.1f}-{ci4k_hi:.1f}%); r_rb={r_rb4k:+.3f} (Z={z4k:+.2f}; pair Z={z4p:+.1f})",
            "test_type": "Wilcoxon signed-rank (khipu-level; asymptotic)",
            "raw_p": p4k,
            "p_display": f"{p4k:.2e}" if p4k > 0 else "p < 1e-300 (asymptotic Z)",
        },
        {
            "regularity": "R5: Vertical descending sibling order (dAP >= 3 cm)",
            "sample_dedup": f"{int(r5_3['n_top_gt_bot'])}/{int(r5_3['n_pairs'])} pairs ({meta5_ded['n_khipus_any']} khipus)",
            "effect_metric": "Frac top>bottom vs dAP=0 neg ctrl",
            "effect_size": f"{100*r5_3['frac_top_gt_bot']:.1f}% (95% CI {100*r5_3['ci_lo']:.1f}-{100*r5_3['ci_hi']:.1f}%) vs {100*r5_0['frac_top_gt_bot']:.1f}% at 0cm (delta={100*(r5_3['frac_top_gt_bot']-r5_0['frac_top_gt_bot']):+.1f} pp)",
            "test_type": "Exact binomial (pair) + khipu Wilcoxon",
            "raw_p": float(r5_3["binom_p"]),
            "p_display": f"{r5_3['binom_p']:.2e} (khipu p={meta5_ded['wilcoxon_p_all']:.2e})",
        },
    ]
    edf = pd.DataFrame(reg_rows)
    # Apply BOTH BH-FDR and Holm-Bonferroni across the 5 candidate regularities!
    pvals = edf["raw_p"].values
    edf["q_bh_5reg"] = false_discovery_control(pvals, method="bh")
    # Holm step-down
    order = np.argsort(pvals)
    holm = np.empty_like(pvals)
    running_max = 0.0
    m_tests = len(pvals)
    for rank, idx in enumerate(order):
        val = min(1.0, pvals[idx] * (m_tests - rank))
        running_max = max(running_max, val)
        holm[idx] = running_max
    edf["p_holm_5reg"] = holm

    # Also save auxiliary R4c exact part-whole & khipu-level summary
    aux_df = pd.DataFrame([
        {"subtest": "R4a_pair_level_dedup", "n": len(g4), "effect": f"frac={100*(g4['cord_val']<g4['p_cord_val']).mean():.1f}%, r_rb={r_rb4p:+.3f}, Z={z4p:+.2f}", "p_bound": f"p < 1e-300 (asymptotic, |Z|={abs(z4p):.1f})" if p4p == 0 else f"{p4p:.2e}"},
        {"subtest": "R4c_part_whole_ge10_dedup", "n": len(m10), "effect": f"obs={obs_pw} vs null_mean={null_pw.mean():.2f} (enrichment={obs_pw/max(null_pw.mean(),1e-9):.2f}x)", "p_bound": f"p={p_pw:.4f} (N_PERM={n_perm_pw}, min=1/{n_perm_pw+1})"},
        {"subtest": "R5_khipu_level_dedup", "n": meta5_ded["n_khipus_any"], "effect": f"gt50={meta5_ded['n_khipus_gt50']}/{meta5_ded['n_khipus_any']} (med={100*meta5_ded['med_frac_all']:.1f}%)", "p_bound": f"Wilcoxon p={meta5_ded['wilcoxon_p_all']:.2e}"},
    ])
    return edf, aux_df


# =====================================================================
# 4. CHECK PART 2 OF UR1175 (GROUPS 22-24, CORDS 105-119) vs LOT VA 63038
# =====================================================================
def check_triptych_part2_vs_va63038() -> pd.DataFrame:
    p192, _ = load_kfg_pendants("UR1175.xlsx")
    part2_vals = np.array([p192[i][0] for i in range(105, 120)], dtype=int)  # 3 groups of 5 = 15 cords
    part3_vals = np.array([p192[i][0] for i in range(120, 225)], dtype=int).reshape(21, 5)
    # Check Ascher Note 5: Part 2 group g (g=0,1,2) vs sum of every 3rd group in Part 3 (7 sets of 3 groups)!
    p3_sums = np.vstack([part3_vals[g::3].sum(axis=0) for g in range(3)])  # (3, 5)
    p2_mat = part2_vals.reshape(3, 5)

    rows = []
    for g in range(3):
        rows.append({
            "check": f"UR1175_Part2_Grp{22+g}_vs_Part3_setsum",
            "target_khipu": "UR1175 (internal Ascher Note 5)",
            "part2_vals": str(p2_mat[g].tolist()),
            "comparison_vals": str(p3_sums[g].tolist()),
            "exact_matches": int((p2_mat[g] == p3_sums[g]).sum()),
            "n_cells": 5,
            "pearson_r": round(float(stats.pearsonr(p2_mat[g], p3_sums[g]).statistic), 4),
            "note": f"Part2 sum={p2_mat[g].sum()} vs Part3 set-sum={p3_sums[g].sum()}",
        })

    # Slide Part 2 (15 values) across ALL 3 khipus in Berlin lot VA 63038: UR231 (KH0468), UR232 (KH0469), UR233 (KH0470)
    for fname, okr, desc in [
        ("UR231.xlsx", "KH0468", "VA 63038(A) Ica/Pisco"),
        ("UR232.xlsx", "KH0469", "VA 63038(B) Huacho"),
        ("UR233.xlsx", "KH0470", "VA 63038 b Huacho"),
    ]:
        pk, _ = load_kfg_pendants(fname)
        arr = np.array([pk[i][0] for i in sorted(pk.keys())], dtype=int)
        # Slide window of 15 (fwd and rev) over raw array AND over non-zero array
        best_ex, best_r, best_pos = -1, -1.0, ""
        for label, seq in [("raw", arr), ("nonzero", arr[arr > 0])]:
            for ori, su in [("fwd", seq), ("rev", seq[::-1])]:
                for st in range(0, max(1, len(su) - 15 + 1)):
                    win = su[st:st + 15]
                    if len(win) < 15:
                        continue
                    ex = int((win == part2_vals).sum())
                    r = float(stats.pearsonr(win, part2_vals).statistic) if np.std(win) > 0 else 0.0
                    if (ex, r) > (best_ex, best_r):
                        best_ex, best_r, best_pos = ex, r, f"{label}_{ori}_pos{st}"
        rows.append({
            "check": f"UR1175_Part2_15cords_vs_{okr}",
            "target_khipu": f"{okr} ({fname[:-5]}, {desc}, n={len(arr)})",
            "part2_vals": str(part2_vals.tolist()),
            "comparison_vals": best_pos,
            "exact_matches": best_ex,
            "n_cells": 15,
            "pearson_r": round(best_r, 4),
            "note": "No Part 2 copy found in VA 63038" if best_ex < 8 else "Potential Part 2 match!",
        })
    return pd.DataFrame(rows)


def main() -> None:
    rng = np.random.default_rng(SEED)
    cohorts = load_cohorts()
    okr_pends = load_okr_pendants()[0] if isinstance(load_okr_pendants(), tuple) else load_okr_pendants()

    print("1/4 Running Step 3 secondary sensitivity (|k_K - k_D| == 1)...", flush=True)
    pm1_df = run_step3_secondary_pm1(cohorts, okr_pends, rng, n_perm=5000)
    pm1_df.to_csv(RESULTS_DIR / "step3_secondary_k_pm1.csv", index=False)
    for coh, g in pm1_df.groupby("cohort"):
        print(f"  {coh}: n_pairs={len(g)}, BH q<0.05={int(g['pass_fdr'].sum())}, max A_JSD={g['A_JSD_pm1'].max():.4f}, min q={g['q_bh'].min():.4f}")

    print("\n2/4 Running Step 3 power analysis (noise curve)...", flush=True)
    pwr_df = run_step3_power_analysis(cohorts, okr_pends, rng, n_trials=120, n_perm_n2=1000)
    pwr_df.to_csv(RESULTS_DIR / "step3_power_analysis.csv", index=False)
    print(pwr_df[pwr_df["noise_sigma"].isin([0.0, 0.10, 0.20, 0.30, 0.50, 1.0])][["k", "noise_sigma", "mean_A_JSD", "power_null1_p05", "power_null2_p05", "power_eval_p05", "power_fdr30_proxy"]].to_string(index=False))

    print("\n3/4 Running deduplicated Phase 2A re-run + effect sizes + 5-test FDR/Holm...", flush=True)
    edf, aux_df = run_phase2a_dedup_audit(rng)
    edf.to_csv(RESULTS_DIR / "phase2a_dedup_effect_sizes.csv", index=False)
    aux_df.to_csv(RESULTS_DIR / "phase2a_dedup_aux_tests.csv", index=False)
    print(edf[["regularity", "sample_dedup", "p_display", "q_bh_5reg", "p_holm_5reg"]].to_string(index=False))
    print("\nAuxiliary dedup tests:")
    print(aux_df.to_string(index=False))

    print("\n4/4 Checking UR1175 Part 2 (Groups 22-24) vs Berlin lot VA 63038...", flush=True)
    p2_df = check_triptych_part2_vs_va63038()
    p2_df.to_csv(RESULTS_DIR / "phase2c_part2_va63038_check.csv", index=False)
    print(p2_df[["check", "target_khipu", "exact_matches", "n_cells", "pearson_r", "note"]].to_string(index=False))


if __name__ == "__main__":
    main()
