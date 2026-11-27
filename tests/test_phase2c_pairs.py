import pandas as pd
from src.data_loader import RESULTS_DIR


def test_berlin_triptych_part1_and_part3_accuracy():
    tdf = pd.read_csv(RESULTS_DIR / "phase2c_berlin_triptych_cells.csv")
    p1 = tdf[tdf["triptych_half"] == "Part1_vs_UR233"]
    p3 = tdf[tdf["triptych_half"] == "Part3_vs_UR232"]
    assert len(p1) == 104
    assert len(p3) == 105
    assert p1["exact_match"].sum() >= 85
    assert p3["exact_match"].sum() >= 68
    # Physical distinctness: master is ~20 cm longer than copy
    assert p1["master_len_cm"].median() - p1["copy_len_cm"].median() > 15.0


def test_reversed_duplicate_kh0121_kh0473_physical_identity():
    dup = pd.read_csv(RESULTS_DIR / "phase2c_reversed_duplicate_kh0121_kh0473.csv")
    assert len(dup) == 29
    assert dup["val_exact"].sum() >= 24
    valid = dup[(dup["kh0121_len_cm"] > 0) & (dup["kh0473_len_cm_rev"] > 0)]
    assert valid["len_abs_diff_cm"].median() <= 2.5


def test_sweep_recovers_positive_controls_and_triptych():
    sw = pd.read_csv(RESULTS_DIR / "phase2c_corpus_pair_sweep.csv")
    cats = set(sw["category"].tolist())
    assert any("Berlin Triptych Part 1" in c for c in cats)
    assert any("Reversed Dup" in c for c in cats)
    assert any("AS125/UR247" in c for c in cats)
    assert any("Puruchuco" in c for c in cats)
