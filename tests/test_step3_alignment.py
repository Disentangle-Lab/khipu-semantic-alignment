"""Unit tests for pre-registered Step 3 colonial alignment functions."""

import numpy as np

from src.step3_colonial_alignment import align_jsd, align_spearman, mean_column_shares, null1_column_perm_p, null2_empirical_p


def test_mean_column_shares_scale_invariance():
    mat1 = np.array([
        [10, 20, 70],
        [20, 40, 140],
    ], dtype=float)
    mat2 = mat1 * 17.5  # exact constant scaling
    p1 = mean_column_shares(mat1)
    p2 = mean_column_shares(mat2)
    assert np.allclose(p1, [0.1, 0.2, 0.7])
    assert np.allclose(p1, p2)


def test_align_jsd_orientation_invariance():
    p_k = np.array([0.5, 0.3, 0.15, 0.05])
    p_d_rev = p_k[::-1]
    score, orient = align_jsd(p_k, p_d_rev)
    assert abs(score - 1.0) < 1e-12
    assert orient == "rev"
    assert abs(align_spearman(p_k, p_d_rev) - 1.0) < 1e-12


def test_null1_exact_resolution_k6():
    # For k=6 with strictly distinct shares, exact match should achieve p = 2 / 720 = 1 / 360
    p_d = np.array([0.40, 0.25, 0.15, 0.10, 0.06, 0.04])
    rng = np.random.default_rng(42)
    p_val = null1_column_perm_p(p_d, p_d, 1.0, rng)
    assert abs(p_val - (2.0 / 720.0)) < 1e-12


def test_null2_empirical_within_row_fallback():
    mat_k = np.tile(np.array([500.0, 200.0, 100.0, 50.0, 10.0, 5.0]), (8, 1))
    p_d = mean_column_shares(mat_k)
    obs, _ = align_jsd(p_d, p_d)
    rng = np.random.default_rng(42)
    # Empty pool triggers within-khipu row-wise column shuffle
    p_val = null2_empirical_p(mat_k, p_d, obs, empirical_pool=[], rng=rng, n_perm=200)
    assert p_val < 0.05
