"""Unit tests for Step 1 and Step 2 core statistical functions."""

import numpy as np
import pandas as pd

from src.step1_archive_invariance import jsd_base2, modal_color_profile, sim_color_modal, sim_numeric_rank
from src.step2_taxonomy_outliers import assign_deterministic_family, compute_x100_groups, subsidiary_slot_concentration


def test_jsd_bounds_and_symmetry():
    p = np.array([0.5, 0.3, 0.2])
    q = np.array([0.1, 0.6, 0.3])
    assert abs(jsd_base2(p, p)) < 1e-12
    assert 0.0 <= jsd_base2(p, q) <= 1.0
    assert abs(jsd_base2(p, q) - jsd_base2(q, p)) < 1e-12
    # Disjoint distributions -> JSD_2 == 1.0
    assert abs(jsd_base2(np.array([1.0, 0.0]), np.array([0.0, 1.0])) - 1.0) < 1e-12


def test_orientation_invariant_similarity():
    c1 = ["W", "B", "R", "G"]
    c2 = ["G", "R", "B", "W"]  # exact reverse
    assert sim_color_modal(c1, c2, orient_inv=True) == 1.0
    assert sim_color_modal(c1, c2, orient_inv=False) == 0.0

    r1 = np.array([1.0, 2.0, 3.0, 4.0])
    r2 = np.array([4.0, 3.0, 2.0, 1.0])
    assert abs(sim_numeric_rank(r1, r2, orient_inv=True) - 1.0) < 1e-12
    assert abs(sim_numeric_rank(r1, r2, orient_inv=False) - (-1.0)) < 1e-12


def test_modal_color_profile_deterministic():
    mat = np.array([
        ["W", "B", "R"],
        ["W", "G", "R"],
        ["B", "B", "__MISSING__"],
    ], dtype=object)
    assert modal_color_profile(mat) == ["W", "B", "R"]


def test_x100_groups_loose_vs_strict():
    mat = np.array([
        [40, 60, 0],     # sum=100, loose only (all mult 10)
        [37, 63, 0],     # sum=100, strict (non-mult 10)
        [100, 0, 0],     # only 1 non-zero -> excluded
        [10, 20, 30],    # sum=60 -> excluded
    ], dtype=float)
    loose, strict = compute_x100_groups(mat)
    assert (loose, strict) == (2, 1)


def test_subsidiary_concentration():
    mat_sub = np.array([
        [5, 0, 0],
        [4, 0, 0],
        [6, 0, 0],
    ])
    tot, ent, pval = subsidiary_slot_concentration(mat_sub)
    assert tot == 15
    assert abs(ent - 0.0) < 1e-12
    assert pval < 1e-4


def test_deterministic_families_and_outliers():
    base = {
        "G": 5, "k": 4, "clean_dec_frac": 1.0, "h5c_sig": False,
        "n_groups_x100_strict": 0, "table_coverage": 0.9,
    }
    assert assign_deterministic_family(pd.Series(base))[0] == "F1_compact_parallel"
    assert assign_deterministic_family(pd.Series({**base, "k": 6}))[0] == "F2_medium_multi_attr"
    assert assign_deterministic_family(pd.Series({**base, "h5c_sig": True}))[0] == "F3_self_totaling"
    assert assign_deterministic_family(pd.Series({**base, "clean_dec_frac": 0.85}))[0] == "F4_nonclean_hybrid"
    assert assign_deterministic_family(pd.Series({**base, "k": 160}))[0] == "OUT"
    assert assign_deterministic_family(pd.Series({**base, "G": 2}))[0] == "OUT"
