"""
step3_colonial_alignment.py — Pre-registered Step 3 (LOCKED to PREREGISTRATION.md SHA-256):
Scale-invariant Table-to-Table Alignment between tabular khipus and structured colonial
Andean khipu transcriptions (`data/colonial_schemas/*.json`).

Strict adherence to frozen PREREGISTRATION.md:
  - Eligibility: exact column-arity match (k_K == k_D)
  - Primary metric: orientation-invariant Jensen-Shannon similarity A_JSD = 1 - min(JSD_2(fwd), JSD_2(rev))
    on row-normalized mean column shares (continuous in [0, 1], weight-free)
  - Secondary robustness metric: orientation-invariant Spearman rank correlation A_rank
  - Resolution-aware dual null (N_PERM = 5,000, BH-FDR q < 0.05):
      * Null 2 (corpus-empirical khipu null): evaluated individually for ALL k >= 3
      * Null 1 (column-order permutation null): evaluated individually via intersection-union
        for k >= 6, and at the aggregate stratum level for k in {3, 4, 5}
  - Primary cohort: n = 86 doubly-tabular khipus (H3a & H3b FDR q < 0.05)
  - Secondary sensitivity cohort (reported strictly separately): H3b-only numeric tabular khipus
  - Cross-check: KFG Excel vertical-position values alongside OKR base-10 values.
"""

from __future__ import annotations

import itertools
import json
import pathlib

import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import false_discovery_control

from src.data_loader import (
    DATA_DIR,
    RESULTS_DIR,
    extract_modal_table,
    load_cohorts,
    load_kfg_excel_pendants,
    load_okr_pendants,
)
from src.step1_archive_invariance import jsd_base2

SCHEMA_DIR = DATA_DIR / "colonial_schemas"
N_PERM = 5000
SEED = 1670


def mean_column_shares(mat: np.ndarray) -> np.ndarray:
    """Normalize each non-zero row to sum to 1, then average across non-zero rows."""
    mat = np.asarray(mat, dtype=float)
    row_sums = mat.sum(axis=1)
    nz = row_sums > 0
    if not np.any(nz):
        return np.full(mat.shape[1], np.nan)
    normed = mat[nz] / row_sums[nz, None]
    return normed.mean(axis=0)


def align_jsd(p_k: np.ndarray, p_d: np.ndarray) -> tuple[float, str]:
    """Primary pre-registered metric A_JSD in [0, 1] and best orientation ('fwd' or 'rev')."""
    d_fwd = jsd_base2(p_k, p_d)
    d_rev = jsd_base2(p_k, p_d[::-1])
    if d_fwd <= d_rev:
        return float(1.0 - d_fwd), "fwd"
    return float(1.0 - d_rev), "rev"


def align_spearman(p_k: np.ndarray, p_d: np.ndarray) -> float:
    """Secondary pre-registered robustness metric A_rank in [-1, 1]."""
    if np.std(p_k) <= 1e-12 or np.std(p_d) <= 1e-12:
        return np.nan
    r_fwd = float(stats.spearmanr(p_k, p_d).statistic)
    r_rev = float(stats.spearmanr(p_k, p_d[::-1]).statistic)
    return max(r_fwd, r_rev)


def jsd_base2_batch(P: np.ndarray, Q: np.ndarray) -> np.ndarray:
    """Vectorized base-2 JSD between rows of P (N, k) and Q (N, k) or broadcastable (1, k)."""
    P = np.asarray(P, dtype=float)
    Q = np.asarray(Q, dtype=float)
    P = P / np.maximum(P.sum(axis=-1, keepdims=True), 1e-15)
    Q = Q / np.maximum(Q.sum(axis=-1, keepdims=True), 1e-15)
    M = 0.5 * (P + Q)
    with np.errstate(divide="ignore", invalid="ignore"):
        kl_p = np.where(P > 0, P * np.log2(P / M), 0.0).sum(axis=-1)
        kl_q = np.where(Q > 0, Q * np.log2(Q / M), 0.0).sum(axis=-1)
    return 0.5 * (kl_p + kl_q)


def align_jsd_batch(P_k: np.ndarray, P_d: np.ndarray) -> np.ndarray:
    """Vectorized orientation-invariant A_JSD = 1 - min(JSD(fwd), JSD(rev))."""
    d_fwd = jsd_base2_batch(P_k, P_d)
    d_rev = jsd_base2_batch(P_k, P_d[..., ::-1])
    return 1.0 - np.minimum(d_fwd, d_rev)


def null1_column_perm_p(
    p_k: np.ndarray, p_d: np.ndarray, obs_jsd: float, rng: np.random.Generator, n_perm: int = N_PERM
) -> float:
    """Exact enumeration for k <= 7; Monte Carlo (n_perm) for k > 7 (vectorized)."""
    k = len(p_d)
    if k <= 7:
        perms = np.array(list(itertools.permutations(range(k))), dtype=int)
        scores = align_jsd_batch(p_k[None, :], p_d[perms])
        return float(np.mean(scores >= obs_jsd - 1e-12))
    idx = np.argsort(rng.random((n_perm, k)), axis=1)
    scores = align_jsd_batch(p_k[None, :], p_d[idx])
    return float((1 + int(np.sum(scores >= obs_jsd - 1e-12))) / (n_perm + 1))


def null2_empirical_p(
    mat_k: np.ndarray,
    p_d: np.ndarray,
    obs_jsd: float,
    empirical_pool: list[np.ndarray],
    rng: np.random.Generator,
    n_perm: int = N_PERM,
) -> float:
    """
    Corpus-empirical khipu null (N_PERM = 5,000, vectorized):
    If >= 20 khipus share k, draw synthetic G x k profiles by sampling G rows from the
    empirical row pool of same-k khipus; otherwise perform within-khipu row-wise column shuffles.
    """
    G, k = mat_k.shape
    if len(empirical_pool) >= 20:
        all_rows = np.vstack(empirical_pool)
        rs = all_rows.sum(axis=1)
        all_norm = all_rows[rs > 0] / rs[rs > 0, None]
        n_norm = all_norm.shape[0]
        # Process in chunks to keep memory modest
        ge = 0
        chunk = 1000
        for start in range(0, n_perm, chunk):
            bs = min(chunk, n_perm - start)
            idx = rng.integers(0, n_norm, size=(bs, G))
            p_syn = all_norm[idx].mean(axis=1)  # (bs, k)
            scores = align_jsd_batch(p_syn, p_d[None, :])
            ge += int(np.sum(scores >= obs_jsd - 1e-12))
        return float((1 + ge) / (n_perm + 1))

    rs = mat_k.sum(axis=1)
    nz_rows = mat_k[rs > 0] / rs[rs > 0, None]
    G_nz = nz_rows.shape[0]
    if G_nz == 0:
        return 1.0
    ge = 0
    chunk = 1000
    for start in range(0, n_perm, chunk):
        bs = min(chunk, n_perm - start)
        p_syn = np.zeros((bs, k), dtype=float)
        for g_i in range(G_nz):
            idx = np.argsort(rng.random((bs, k)), axis=1)
            p_syn += nz_rows[g_i][idx]
        p_syn /= G_nz
        scores = align_jsd_batch(p_syn, p_d[None, :])
        ge += int(np.sum(scores >= obs_jsd - 1e-12))
    return float((1 + ge) / (n_perm + 1))


def load_colonial_schemas(schema_dir: pathlib.Path = SCHEMA_DIR) -> list[dict]:
    schemas = []
    for path in sorted(schema_dir.glob("*.json")):
        with open(path, encoding="utf-8") as f:
            obj = json.load(f)
        mat = np.asarray(obj["value_matrix"], dtype=float)
        obj["_p_d"] = mean_column_shares(mat)
        obj["_k"] = mat.shape[1]
        obj["_m"] = mat.shape[0]
        schemas.append(obj)
    return schemas


def build_empirical_pools(pendants: pd.DataFrame) -> dict[int, list[np.ndarray]]:
    pools: dict[int, list[np.ndarray]] = {}
    for _, g in pendants.groupby("KHIPU_ID", sort=False):
        if len(g) < 10:
            continue
        tab = extract_modal_table(g)
        if tab is not None:
            pools.setdefault(tab["k"], []).append(tab["mat_val"])
    return pools


def run_cohort_alignment(
    cohort_df: pd.DataFrame,
    cohort_name: str,
    pend_by_kid: dict[int, pd.DataFrame],
    schemas: list[dict],
    emp_pools: dict[int, list[np.ndarray]],
    rng: np.random.Generator,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows = []
    agg_rows = []

    for sc in schemas:
        k_d = sc["_k"]
        p_d = sc["_p_d"]
        pool_k = emp_pools.get(k_d, [])
        pool_mode = f"empirical_corpus(n={len(pool_k)})" if len(pool_k) >= 20 else f"within_row_shuf(n_pool={len(pool_k)})"

        elig_pks = []
        elig_obs = []

        for r in cohort_df.itertuples(index=False):
            if int(r.h3_m) != k_d:
                continue
            kid = int(r.KHIPU_ID)
            tab = extract_modal_table(pend_by_kid[kid])
            if tab is None or tab["k"] != k_d:
                continue
            mat_k = tab["mat_val"]
            p_k = mean_column_shares(mat_k)
            if not np.all(np.isfinite(p_k)):
                continue

            obs_jsd, orient = align_jsd(p_k, p_d)
            obs_rho = align_spearman(p_k, p_d)
            p1 = null1_column_perm_p(p_k, p_d, obs_jsd, rng)
            p2 = null2_empirical_p(mat_k, p_d, obs_jsd, pool_k, rng)

            # Independent KFG Excel value cross-check when exact KFG structure matches
            kfg_jsd = np.nan
            lf = r.legal_filename if isinstance(r.legal_filename, str) else ""
            if lf:
                kp = load_kfg_excel_pendants(lf)
                if kp is not None and len(kp) > 0:
                    kp_cl = kp[kp["kfg_cluster"] >= 0]
                    csz = kp_cl.groupby("kfg_cluster", sort=False).size()
                    valid = csz[csz == k_d]
                    if len(valid) >= 2:
                        kmat = np.array([kp_cl[kp_cl["kfg_cluster"] == c]["kfg_val"].values.astype(float) for c in valid.index])
                        pk_kfg = mean_column_shares(kmat)
                        if np.all(np.isfinite(pk_kfg)):
                            kfg_jsd = align_jsd(pk_kfg, p_d)[0]

            elig_pks.append(p_k)
            elig_obs.append(obs_jsd)
            rows.append({
                "cohort": cohort_name,
                "schema_id": sc["schema_id"],
                "k": k_d,
                "OKR_NUM": r.OKR_NUM,
                "kfg_name": r.kfg_name if pd.notna(r.kfg_name) else "",
                "provenance": r.provenance_clean,
                "cataloger": r.cataloger,
                "G": tab["G"],
                "clean_dec_frac": round(float(r.clean_dec_frac), 3),
                "A_JSD": round(obs_jsd, 4),
                "A_JSD_kfg": round(kfg_jsd, 4) if np.isfinite(kfg_jsd) else np.nan,
                "orient": orient,
                "A_rank": round(obs_rho, 4) if np.isfinite(obs_rho) else np.nan,
                "p_null1": round(p1, 5),
                "p_null2": round(p2, 5),
                "p_both": round(max(p1, p2), 5),
                "null2_mode": pool_mode,
                "khipu_shares": np.round(p_k, 3).tolist(),
                "schema_shares": np.round(p_d, 3).tolist(),
            })

        # Aggregate Null 1 stratum test (vectorized; mandatory per PREREGISTRATION for k <= 5; reported for all k)
        if elig_pks:
            obs_mean = float(np.mean(elig_obs))
            P_elig = np.array(elig_pks, dtype=float)  # (n_elig, k)
            if k_d <= 7:
                perms = np.array(list(itertools.permutations(range(k_d))), dtype=int)
                D_perm = p_d[perms]  # (P, k)
                # Compute (P, n_elig) scores via broadcasting
                sc_mat = align_jsd_batch(P_elig[None, :, :], D_perm[:, None, :])
                null_means = sc_mat.mean(axis=1)
                p_agg1 = float(np.mean(null_means >= obs_mean - 1e-12))
            else:
                idx = np.argsort(rng.random((N_PERM, k_d)), axis=1)
                D_perm = p_d[idx]
                sc_mat = align_jsd_batch(P_elig[None, :, :], D_perm[:, None, :])
                null_means = sc_mat.mean(axis=1)
                p_agg1 = float((1 + int(np.sum(null_means >= obs_mean - 1e-12))) / (N_PERM + 1))
            agg_rows.append({
                "cohort": cohort_name,
                "schema_id": sc["schema_id"],
                "k": k_d,
                "n_khipus_tested": len(elig_pks),
                "mean_A_JSD": round(obs_mean, 4),
                "max_A_JSD": round(float(np.max(elig_obs)), 4),
                "p_agg_null1": round(p_agg1, 5),
                "null2_mode": pool_mode,
            })

    df_pairs = pd.DataFrame(rows)
    if len(df_pairs) > 0:
        df_pairs["q_null2"] = false_discovery_control(df_pairs["p_null2"].values, method="bh")
        df_pairs["q_both_k6plus"] = np.nan
        m6 = df_pairs["k"] >= 6
        if m6.sum() > 0:
            df_pairs.loc[m6, "q_both_k6plus"] = false_discovery_control(df_pairs.loc[m6, "p_both"].values, method="bh")
        df_pairs["passes_prereg"] = np.where(
            df_pairs["k"] >= 6,
            df_pairs["q_both_k6plus"] < 0.05,
            df_pairs["q_null2"] < 0.05,
        )
    return df_pairs, pd.DataFrame(agg_rows)


def main():
    schemas = load_colonial_schemas()
    print(f"Loaded {len(schemas)} verified colonial schemas:")
    for sc in schemas:
        print(f"  {sc['schema_id']:28s} ({sc['_m']:2d}x{sc['_k']:2d}) | shares={np.round(sc['_p_d'], 3).tolist()}")

    # Text-to-text positive control: Huaraz 1532 (m_60) vs Huaraz 1558 (m_63) (both k=7!)
    sc_map = {s["schema_id"]: s for s in schemas}
    if "huaraz_1532_full_k7" in sc_map and "huaraz_1558_full_k7" in sc_map:
        p32 = sc_map["huaraz_1532_full_k7"]["_p_d"]
        p58 = sc_map["huaraz_1558_full_k7"]["_p_d"]
        a_jsd, ori = align_jsd(p32, p58)
        a_rho = align_spearman(p32, p58)
        p1 = null1_column_perm_p(p32, p58, a_jsd, np.random.default_rng(SEED))
        print(
            f"\n=== TEXT-TO-TEXT CALIBRATION: Huaraz 1532 (m_60) vs Huaraz 1558 (m_63) (k=7) ===\n"
            f"  A_JSD = {a_jsd:.4f} ({ori}) | A_rank = {a_rho:.4f} | exact Null 1 p = {p1:.5f} (floor = 2/5040 = 0.00040)"
        )

    print("\nLoading OKR v2.1.0 pendants, empirical pools, and cohorts...")
    pendants, _ = load_okr_pendants()
    cohorts = load_cohorts()
    pend_by_kid = {int(k): g for k, g in pendants.groupby("KHIPU_ID", sort=False)}
    emp_pools = build_empirical_pools(pendants)

    dbl = cohorts[cohorts["doubly_tabular"]].copy()
    h3b_only = cohorts[cohorts["h3b_only"]].copy()
    print(f"Primary cohort: doubly_tabular={len(dbl)} | Secondary sensitivity: H3b-only={len(h3b_only)}")

    rng = np.random.default_rng(SEED)
    p_prim, a_prim = run_cohort_alignment(dbl, "primary_86", pend_by_kid, schemas, emp_pools, rng)
    p_sens, a_sens = run_cohort_alignment(h3b_only, "sens_H3b_only", pend_by_kid, schemas, emp_pools, rng)

    all_pairs = pd.concat([p_prim, p_sens], ignore_index=True)
    all_agg = pd.concat([a_prim, a_sens], ignore_index=True)
    all_pairs.to_csv(RESULTS_DIR / "step3_alignment_pairs.csv", index=False)
    all_agg.to_csv(RESULTS_DIR / "step3_alignment_stratum_agg.csv", index=False)

    print("\n=== STEP 3 STRATUM AGGREGATE SUMMARY (PRIMARY 86 & SENSITIVITY H3b-ONLY) ===")
    print(all_agg.to_string(index=False))

    print("\n=== STEP 3 PRIMARY 86: PAIRS PASSING PRE-REGISTERED FDR (q < 0.05) ===")
    hits_p = p_prim[p_prim["passes_prereg"]] if len(p_prim) else pd.DataFrame()
    print(f"Total tested pairs in primary_86: {len(p_prim)} | Passing pre-reg FDR: {len(hits_p)}")
    if len(hits_p):
        print(
            hits_p.sort_values("A_JSD", ascending=False)[
                ["schema_id", "k", "OKR_NUM", "kfg_name", "provenance", "G", "A_JSD", "A_JSD_kfg", "orient", "A_rank", "p_null1", "p_null2", "q_null2", "q_both_k6plus"]
            ].round(4).to_string(index=False)
        )
    else:
        print("  (0 pairs passed pre-registered FDR in primary_86 — reporting null result transparently!)")

    print("\n=== STEP 3 PRIMARY 86: TOP 10 PAIRS BY A_JSD (DIAGNOSTIC VIEW ONLY) ===")
    if len(p_prim):
        print(
            p_prim.sort_values("A_JSD", ascending=False).head(10)[
                ["schema_id", "k", "OKR_NUM", "kfg_name", "provenance", "G", "A_JSD", "A_JSD_kfg", "orient", "A_rank", "p_null1", "p_null2", "q_null2", "q_both_k6plus"]
            ].round(4).to_string(index=False)
        )

    print("\n=== STEP 3 SECONDARY SENSITIVITY (H3b-ONLY): PAIRS PASSING PRE-REG FDR ===")
    hits_s = p_sens[p_sens["passes_prereg"]] if len(p_sens) else pd.DataFrame()
    print(f"Total tested pairs in sens_H3b_only: {len(p_sens)} | Passing pre-reg FDR: {len(hits_s)}")
    if len(hits_s):
        print(
            hits_s.sort_values("A_JSD", ascending=False)[
                ["schema_id", "k", "OKR_NUM", "kfg_name", "provenance", "G", "A_JSD", "A_JSD_kfg", "orient", "A_rank", "p_null1", "p_null2", "q_null2", "q_both_k6plus"]
            ].round(4).to_string(index=False)
        )


if __name__ == "__main__":
    main()
