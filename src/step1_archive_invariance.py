"""
step1_archive_invariance.py — Pre-registered Step 1:
Inter-khipu archive invariance of column-color (S1a) and column-numeric rank (S1b) profiles,
strictly controlled for column count (k = h3_m), archive sample size, and cataloger.

Implements PREREGISTRATION.md Section 1:
  - Orientation-invariant pairwise similarity within exact k (forward-only in supplement)
  - Macro-averaged excess T_macro (equal weight 1/|S| per (site, k) cell with n >= 2)
    alongside micro-averaged (pair-weighted) excess
  - Strict k-stratified site-label permutation (N_PERM = 5,000)
  - Same-(k, cataloger) stratified permutation control
  - Matched-size subsampling (n_sub = 2, 1,000 draws) for every cell with n > 2
  - Leave-one-archive-out sensitivity
  - Robustness: full column color distribution JSD similarity & deduplicated / same-object collapsed corpus.
"""

from __future__ import annotations

import itertools
from collections import Counter

import numpy as np
import pandas as pd
from scipy import stats

from src.data_loader import (
    RESULTS_DIR,
    extract_modal_table,
    load_cohorts,
    load_kfg_excel_pendants,
    load_okr_pendants,
)

N_PERM = 5000
N_SUB_DRAWS = 1000
SEED = 2026


def jsd_base2(p: np.ndarray, q: np.ndarray) -> float:
    p = np.asarray(p, dtype=float)
    q = np.asarray(q, dtype=float)
    sp, sq = p.sum(), q.sum()
    if sp <= 0 or sq <= 0:
        return np.nan
    p = p / sp
    q = q / sq
    m = 0.5 * (p + q)
    def kl(a, b):
        mask = a > 0
        return float(np.sum(a[mask] * np.log2(a[mask] / b[mask])))
    return 0.5 * kl(p, m) + 0.5 * kl(q, m)


def modal_color_profile(mat_col: np.ndarray) -> list[str]:
    """Deterministic modal FULL_COLOR per column slot 1..k (excluding __MISSING__ if possible)."""
    _, k = mat_col.shape
    prof = []
    for j in range(k):
        col = [c for c in mat_col[:, j] if c != "__MISSING__"]
        if not col:
            prof.append("__MISSING__")
        else:
            cnt = Counter(col)
            # Sort by (-count, color_str) for deterministic tie-breaking
            best = sorted(cnt.items(), key=lambda x: (-x[1], x[0]))[0][0]
            prof.append(best)
    return prof


def col_color_dists(mat_col: np.ndarray) -> list[dict[str, float]]:
    """Per-column color probability distribution for JSD robustness check."""
    _, k = mat_col.shape
    dists = []
    for j in range(k):
        col = [c for c in mat_col[:, j] if c != "__MISSING__"]
        if not col:
            dists.append({})
        else:
            cnt = Counter(col)
            tot = float(len(col))
            dists.append({c: v / tot for c, v in cnt.items()})
    return dists


def sim_color_modal(c1: list[str], c2: list[str], orient_inv: bool = True) -> float:
    k = len(c1)
    def match_rate(a, b):
        valid = [(x, y) for x, y in zip(a, b) if x != "__MISSING__" and y != "__MISSING__"]
        if not valid:
            return np.nan
        return sum(x == y for x, y in valid) / len(valid)
    fwd = match_rate(c1, c2)
    if not orient_inv:
        return fwd
    rev = match_rate(c1, c2[::-1])
    vals = [v for v in (fwd, rev) if np.isfinite(v)]
    return max(vals) if vals else np.nan


def sim_color_jsd(d1: list[dict[str, float]], d2: list[dict[str, float]], orient_inv: bool = True) -> float:
    def score(a, b):
        slot_sims = []
        for da, db in zip(a, b):
            if not da or not db:
                continue
            keys = sorted(set(da) | set(db))
            pa = np.array([da.get(k, 0.0) for k in keys])
            pb = np.array([db.get(k, 0.0) for k in keys])
            slot_sims.append(1.0 - jsd_base2(pa, pb))
        return float(np.mean(slot_sims)) if slot_sims else np.nan
    fwd = score(d1, d2)
    if not orient_inv:
        return fwd
    rev = score(d1, d2[::-1])
    vals = [v for v in (fwd, rev) if np.isfinite(v)]
    return max(vals) if vals else np.nan


def sim_numeric_rank(r1: np.ndarray, r2: np.ndarray, orient_inv: bool = True) -> float:
    def corr(a, b):
        if np.std(a) <= 1e-12 or np.std(b) <= 1e-12:
            return np.nan
        return float(stats.spearmanr(a, b).statistic)
    fwd = corr(r1, r2)
    if not orient_inv:
        return fwd
    rev = corr(r1, r2[::-1])
    vals = [v for v in (fwd, rev) if np.isfinite(v)]
    return max(vals) if vals else np.nan


def build_profiles(pendants: pd.DataFrame, cohorts: pd.DataFrame, target_kids: set[int] | None = None) -> dict[int, dict]:
    """Extract column-color and column-numeric profiles for target tabular khipus."""
    profiles = {}
    legal_map = cohorts.set_index("KHIPU_ID")["legal_filename"].to_dict()

    for kid, g in pendants.groupby("KHIPU_ID", sort=False):
        if target_kids is not None and int(kid) not in target_kids:
            continue
        if len(g) < 10:
            continue
        tab = extract_modal_table(g)
        if tab is None:
            continue
        log_mat = np.log1p(tab["mat_val"])
        col_means = log_mat.mean(axis=0)
        # Also compute KFG vertical-value column means when exact KFG cluster structure matches (G, k)
        kfg_col_means = None
        lf = legal_map.get(kid)
        if isinstance(lf, str) and lf:
            kp = load_kfg_excel_pendants(lf)
            if kp is not None and len(kp) > 0:
                kp_cl = kp[kp["kfg_cluster"] >= 0]
                csz = kp_cl.groupby("kfg_cluster", sort=False).size()
                valid = csz[csz == tab["k"]]
                if len(valid) >= 2:
                    cids = valid.index.tolist()
                    kmat = np.array([kp_cl[kp_cl["kfg_cluster"] == c]["kfg_val"].values.astype(float) for c in cids])
                    kfg_col_means = np.log1p(kmat).mean(axis=0)

        profiles[int(kid)] = {
            "G": tab["G"],
            "k": tab["k"],
            "modal_colors": modal_color_profile(tab["mat_col"]),
            "color_dists": col_color_dists(tab["mat_col"]),
            "col_means": col_means,
            "kfg_col_means": kfg_col_means,
        }
    return profiles


def run_archive_test(
    df_sub: pd.DataFrame,
    profiles: dict[int, dict],
    metric_name: str,
    sim_fn,
    stratify_cataloger: bool = False,
    n_perm: int = N_PERM,
    seed: int = SEED,
) -> tuple[dict, pd.DataFrame, pd.DataFrame]:
    """
    Run strict k-stratified (or (k, cataloger)-stratified) site-label permutation test.
    Computes macro-averaged T_macro (equal weight per (site, k) cell with n>=2) AND micro-averaged excess.
    """
    rng = np.random.default_rng(seed)

    # Filter to known provenance & available profile
    rows = []
    for r in df_sub.itertuples(index=False):
        kid = int(r.KHIPU_ID)
        if kid not in profiles or r.provenance_clean == "Unknown":
            continue
        rows.append((kid, r.OKR_NUM, r.provenance_clean, r.cataloger, profiles[kid]["k"]))

    if not rows:
        return {"metric": metric_name, "n_khipus": 0}, pd.DataFrame(), pd.DataFrame()

    meta = pd.DataFrame(rows, columns=["KHIPU_ID", "OKR_NUM", "site", "cataloger", "k"])

    # Keep only k strata that have at least one site with n>=2 AND at least 2 distinct sites (so between-site exists)
    valid_ks = []
    for k_val, grp in meta.groupby("k"):
        sc = grp["site"].value_counts()
        if len(sc) >= 2 and sc.iloc[0] >= 2:
            valid_ks.append(k_val)
    meta = meta[meta["k"].isin(valid_ks)].reset_index(drop=True)
    if len(meta) == 0:
        return {"metric": metric_name, "n_khipus": 0}, pd.DataFrame(), pd.DataFrame()

    # Precompute pairwise similarity matrix per k stratum
    stratum_data = {}
    pair_records = []
    for k_val, grp in meta.groupby("k"):
        grp = grp.reset_index(drop=True)
        n_k = len(grp)
        S = np.full((n_k, n_k), np.nan, dtype=float)
        kids = grp["KHIPU_ID"].values
        okrs = grp["OKR_NUM"].values
        sites = grp["site"].values
        cats = grp["cataloger"].values
        for i in range(n_k):
            for j in range(i + 1, n_k):
                val = sim_fn(profiles[int(kids[i])], profiles[int(kids[j])])
                S[i, j] = S[j, i] = val
                pair_records.append({
                    "metric": metric_name,
                    "k": int(k_val),
                    "khipu_a": okrs[i],
                    "khipu_b": okrs[j],
                    "site_a": sites[i],
                    "site_b": sites[j],
                    "same_site": bool(sites[i] == sites[j]),
                    "same_cataloger": bool(cats[i] == cats[j]),
                    "similarity": val,
                })
        iu = np.triu_indices(n_k, k=1)
        uniq_sites, site_codes = np.unique(sites, return_inverse=True)
        _, cat_codes = np.unique(cats, return_inverse=True)
        cat_groups = [np.where(cat_codes == c)[0] for c in range(int(cat_codes.max()) + 1)]
        s_vec = S[iu]
        valid = np.isfinite(s_vec)
        stratum_data[int(k_val)] = {
            "sites": sites,
            "uniq_sites": uniq_sites,
            "site_codes": site_codes.astype(np.int32),
            "n_sites": len(uniq_sites),
            "cat_groups": [g for g in cat_groups if len(g) > 1],
            "iu0": iu[0][valid],
            "iu1": iu[1][valid],
            "s_vec": s_vec[valid],
            "n": n_k,
        }

    def eval_stats_fast(code_map: dict[int, np.ndarray], collect_cells: bool = False):
        cell_sum = 0.0
        cell_cnt = 0
        w_sum = 0.0
        w_cnt = 0
        b_sum = 0.0
        b_cnt = 0
        cell_rows = []
        for k_val, sd in stratum_data.items():
            sc = code_map[k_val]
            sa = sc[sd["iu0"]]
            sb = sc[sd["iu1"]]
            sv = sd["s_vec"]
            if len(sv) == 0:
                continue
            same = sa == sb
            diff = ~same
            n_diff = int(diff.sum())
            if n_diff == 0:
                continue
            sv_diff = sv[diff]
            b_k_sum = float(sv_diff.sum())
            mean_between_k = b_k_sum / n_diff
            b_sum += b_k_sum
            b_cnt += n_diff

            n_same = int(same.sum())
            if n_same == 0:
                continue
            sa_same = sa[same]
            sv_same = sv[same]
            w_sum += float(sv_same.sum())
            w_cnt += n_same

            n_s = sd["n_sites"]
            cnt_s = np.bincount(sa_same, minlength=n_s)
            sum_s = np.bincount(sa_same, weights=sv_same, minlength=n_s)
            has_s = cnt_s > 0
            means_s = sum_s[has_s] / cnt_s[has_s]
            cell_sum += float(np.sum(means_s - mean_between_k))
            cell_cnt += int(has_s.sum())

            if collect_cells:
                uniq_s = sd["uniq_sites"]
                for s_idx in np.where(has_s)[0]:
                    w_cell = float(sum_s[s_idx] / cnt_s[s_idx])
                    cell_rows.append({
                        "metric": metric_name,
                        "site": str(uniq_s[s_idx]),
                        "k": k_val,
                        "n_khipus_cell": int(np.sum(sc == s_idx)),
                        "n_pairs_within": int(cnt_s[s_idx]),
                        "mean_within": w_cell,
                        "mean_between_k": mean_between_k,
                        "excess": w_cell - mean_between_k,
                    })
        t_macro = (cell_sum / cell_cnt) if cell_cnt > 0 else np.nan
        t_micro = ((w_sum / w_cnt) - (b_sum / b_cnt)) if (w_cnt > 0 and b_cnt > 0) else np.nan
        return t_macro, t_micro, cell_rows

    obs_map = {k_val: sd["site_codes"] for k_val, sd in stratum_data.items()}
    obs_macro, obs_micro, obs_cells = eval_stats_fast(obs_map, collect_cells=True)

    null_macro = np.empty(n_perm, dtype=float)
    null_micro = np.empty(n_perm, dtype=float)
    for b in range(n_perm):
        perm_map = {}
        for k_val, sd in stratum_data.items():
            sc = sd["site_codes"]
            if not stratify_cataloger:
                perm_map[k_val] = rng.permutation(sc)
            else:
                out_sc = sc.copy()
                for ix in sd["cat_groups"]:
                    out_sc[ix] = rng.permutation(sc[ix])
                perm_map[k_val] = out_sc
        tm, tmi, _ = eval_stats_fast(perm_map, collect_cells=False)
        null_macro[b] = tm
        null_micro[b] = tmi

    valid_nm = null_macro[np.isfinite(null_macro)]
    valid_nmi = null_micro[np.isfinite(null_micro)]
    p_macro = float((1 + np.sum(valid_nm >= obs_macro)) / (len(valid_nm) + 1)) if len(valid_nm) else np.nan
    p_micro = float((1 + np.sum(valid_nmi >= obs_micro)) / (len(valid_nmi) + 1)) if len(valid_nmi) else np.nan

    summary = {
        "metric": metric_name,
        "stratify_cataloger": stratify_cataloger,
        "n_khipus": int(len(meta)),
        "n_k_strata": int(len(stratum_data)),
        "n_cells": int(len(obs_cells)),
        "T_macro_obs": round(obs_macro, 4) if np.isfinite(obs_macro) else np.nan,
        "T_macro_null_mean": round(float(np.mean(valid_nm)), 4) if len(valid_nm) else np.nan,
        "p_macro": round(p_macro, 4) if np.isfinite(p_macro) else np.nan,
        "T_micro_obs": round(obs_micro, 4) if np.isfinite(obs_micro) else np.nan,
        "T_micro_null_mean": round(float(np.mean(valid_nmi)), 4) if len(valid_nmi) else np.nan,
        "p_micro": round(p_micro, 4) if np.isfinite(p_micro) else np.nan,
    }
    return summary, pd.DataFrame(obs_cells), pd.DataFrame(pair_records)


def run_matched_subsampling(
    df_sub: pd.DataFrame,
    profiles: dict[int, dict],
    metric_name: str,
    sim_fn,
    n_draws: int = N_SUB_DRAWS,
    seed: int = SEED + 1,
) -> pd.DataFrame:
    """For every (site, k) cell with n > 2, subsample n_sub = 2 khipus repeatedly (1,000 draws)."""
    rng = np.random.default_rng(seed)
    rows = []
    for r in df_sub.itertuples(index=False):
        kid = int(r.KHIPU_ID)
        if kid in profiles and r.provenance_clean != "Unknown":
            rows.append((kid, r.provenance_clean, profiles[kid]["k"]))
    meta = pd.DataFrame(rows, columns=["KHIPU_ID", "site", "k"])
    out = []
    for (site, k_val), grp in meta.groupby(["site", "k"]):
        kids = grp["KHIPU_ID"].values
        if len(kids) < 2:
            continue
        # Between-site baseline in same k
        others = meta[(meta["k"] == k_val) & (meta["site"] != site)]["KHIPU_ID"].values
        if len(others) == 0:
            continue
        bet = [
            sim_fn(profiles[int(a)], profiles[int(b)])
            for a in kids for b in others
        ]
        bet = [v for v in bet if np.isfinite(v)]
        if not bet:
            continue
        mean_bet = float(np.mean(bet))

        # Full cell within-site
        all_pairs = [
            sim_fn(profiles[int(a)], profiles[int(b)])
            for a, b in itertools.combinations(kids, 2)
        ]
        all_pairs = [v for v in all_pairs if np.isfinite(v)]
        if not all_pairs:
            continue
        full_mean = float(np.mean(all_pairs))

        # Matched n_sub = 2 draws (each draw picks 2 distinct khipus -> 1 pair)
        if len(kids) > 2:
            draws = rng.choice(all_pairs, size=n_draws, replace=True)
            lo, med, hi = np.percentile(draws, [2.5, 50.0, 97.5])
            frac_gt_bet = float(np.mean(draws > mean_bet))
        else:
            lo = med = hi = full_mean
            frac_gt_bet = float(full_mean > mean_bet)

        out.append({
            "metric": metric_name,
            "site": site,
            "k": int(k_val),
            "n_cell": int(len(kids)),
            "n_pairs": int(len(all_pairs)),
            "full_within_mean": round(full_mean, 4),
            "between_k_mean": round(mean_bet, 4),
            "full_excess": round(full_mean - mean_bet, 4),
            "sub2_median": round(float(med), 4),
            "sub2_ci_lo": round(float(lo), 4),
            "sub2_ci_hi": round(float(hi), 4),
            "sub2_frac_gt_between": round(frac_gt_bet, 3),
        })
    return pd.DataFrame(out)


def main():
    print("Loading OKR v2.1.0 pendants and cohorts...")
    pendants, _ = load_okr_pendants()
    cohorts = load_cohorts()
    dbl = cohorts[cohorts["doubly_tabular"]].copy()
    dbl_dedup = dbl[dbl["keep_dedup"]].copy()
    h3a = cohorts[cohorts["h3a_sig"]].copy()
    h3b = cohorts[cohorts["h3b_sig"]].copy()
    target_kids = set(dbl["KHIPU_ID"].astype(int)) | set(h3a["KHIPU_ID"].astype(int)) | set(h3b["KHIPU_ID"].astype(int))
    profiles = build_profiles(pendants, cohorts, target_kids=target_kids)

    print(
        f"Cohorts ready: doubly_tabular={len(dbl)} (dedup={len(dbl_dedup)}), "
        f"H3a_sig={len(h3a)}, H3b_sig={len(h3b)}"
    )

    sim_specs = [
        ("S1a_color_modal_orient_inv", lambda a, b: sim_color_modal(a["modal_colors"], b["modal_colors"], True)),
        ("S1a_color_modal_fwd_only",   lambda a, b: sim_color_modal(a["modal_colors"], b["modal_colors"], False)),
        ("S1a_color_jsd_orient_inv",   lambda a, b: sim_color_jsd(a["color_dists"], b["color_dists"], True)),
        ("S1b_num_rank_okr_orient_inv", lambda a, b: sim_numeric_rank(a["col_means"], b["col_means"], True)),
        ("S1b_num_rank_okr_fwd_only",  lambda a, b: sim_numeric_rank(a["col_means"], b["col_means"], False)),
        (
            "S1b_num_rank_kfg_orient_inv",
            lambda a, b: sim_numeric_rank(a["kfg_col_means"], b["kfg_col_means"], True)
            if (a["kfg_col_means"] is not None and b["kfg_col_means"] is not None)
            else np.nan,
        ),
    ]

    all_tests = []
    all_cells = []
    all_pairs = []
    all_sub = []

    # 1. Primary cohort (n=86) — raw, dedup, and cataloger-stratified
    for tag, sub_df, strat_cat in [
        ("primary_86", dbl, False),
        ("primary_86_dedup", dbl_dedup, False),
        ("primary_86_cat_strat", dbl, True),
    ]:
        for mname, fn in sim_specs:
            # Run cataloger-stratified & dedup on primary metrics; full spec on primary_86
            if tag != "primary_86" and "fwd_only" in mname:
                continue
            s, cells, pairs = run_archive_test(sub_df, profiles, f"{tag}__{mname}", fn, stratify_cataloger=strat_cat)
            s["cohort"] = tag
            all_tests.append(s)
            if tag == "primary_86":
                all_cells.append(cells)
                if mname in ("S1a_color_modal_orient_inv", "S1b_num_rank_okr_orient_inv"):
                    all_pairs.append(pairs)
                    all_sub.append(run_matched_subsampling(sub_df, profiles, mname, fn))

    # 2. Leave-one-archive-out sensitivity on primary_86 (Incahuasi, Pachacamac)
    for drop_site in ("Incahuasi", "Pachacamac"):
        sub_loo = dbl[dbl["provenance_clean"] != drop_site].copy()
        for mname, fn in [
            ("S1a_color_modal_orient_inv", sim_specs[0][1]),
            ("S1a_color_jsd_orient_inv", sim_specs[2][1]),
            ("S1b_num_rank_okr_orient_inv", sim_specs[3][1]),
        ]:
            s, _, _ = run_archive_test(sub_loo, profiles, f"LOO_no_{drop_site}__{mname}", fn)
            s["cohort"] = f"LOO_no_{drop_site}"
            all_tests.append(s)

    # 3. Expanded sensitivity cohorts (H3a_sig n=133 for color; H3b_sig n=133 for numeric)
    for tag, sub_df, mname, fn in [
        ("sens_H3a_133", h3a, "S1a_color_modal_orient_inv", sim_specs[0][1]),
        ("sens_H3a_133", h3a, "S1a_color_jsd_orient_inv", sim_specs[2][1]),
        ("sens_H3a_133_dedup", h3a[h3a["keep_dedup"]], "S1a_color_modal_orient_inv", sim_specs[0][1]),
        ("sens_H3b_133", h3b, "S1b_num_rank_okr_orient_inv", sim_specs[3][1]),
        ("sens_H3b_133_dedup", h3b[h3b["keep_dedup"]], "S1b_num_rank_okr_orient_inv", sim_specs[3][1]),
    ]:
        s, cells, _ = run_archive_test(sub_df, profiles, f"{tag}__{mname}", fn)
        s["cohort"] = tag
        all_tests.append(s)
        if tag in ("sens_H3a_133", "sens_H3b_133"):
            all_cells.append(cells)

    tests_df = pd.DataFrame(all_tests)
    cells_df = pd.concat(all_cells, ignore_index=True) if all_cells else pd.DataFrame()
    pairs_df = pd.concat(all_pairs, ignore_index=True) if all_pairs else pd.DataFrame()
    sub_df_out = pd.concat(all_sub, ignore_index=True) if all_sub else pd.DataFrame()

    tests_df.to_csv(RESULTS_DIR / "step1_archive_tests.csv", index=False)
    cells_df.to_csv(RESULTS_DIR / "step1_cell_summary.csv", index=False)
    pairs_df.to_csv(RESULTS_DIR / "step1_pair_details.csv", index=False)
    sub_df_out.to_csv(RESULTS_DIR / "step1_matched_subsample.csv", index=False)

    print("\n=== STEP 1 SUMMARY: INTER-KHIPU ARCHIVE INVARIANCE (N_PERM=5,000) ===")
    cols_show = ["metric", "n_khipus", "n_cells", "T_macro_obs", "p_macro", "T_micro_obs", "p_micro"]
    print(tests_df[cols_show].to_string(index=False))

    print("\n=== STEP 1 CELL BREAKDOWN (PRIMARY 86: S1a Color Modal & S1b Numeric Rank) ===")
    c_prim = cells_df[cells_df["metric"].isin([
        "primary_86__S1a_color_modal_orient_inv",
        "primary_86__S1b_num_rank_okr_orient_inv",
    ])]
    print(c_prim.round(4).to_string(index=False))

    print("\n=== STEP 1 MATCHED-SIZE SUBSAMPLING (n_sub=2, 1,000 draws) ===")
    print(sub_df_out.to_string(index=False))


if __name__ == "__main__":
    main()
