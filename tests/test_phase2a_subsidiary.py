import sqlite3
import pandas as pd
import pytest
from src.data_loader import find_okr_db, KFG_XLSX_DIR, RESULTS_DIR, UPSTREAM_DIR
from src.phase2a_subsidiary_grammar import canon_strict, get_disjoint_kfg_files


@pytest.mark.requires_raw_data
def test_kfg_oos_strictly_disjoint_from_okr():
    try:
        db_path = find_okr_db()
    except FileNotFoundError:
        pytest.skip("Local OKR database not present (run scripts/fetch_data.sh)")
    if not KFG_XLSX_DIR.exists():
        pytest.skip("Local KFG Excel workbooks not present (run scripts/fetch_data.sh)")
    con = sqlite3.connect(db_path)
    okr_m = pd.read_sql("SELECT KHIPU_ID, OKR_NUM, INVESTIGATOR_NUM, orig_inv_num, MUSEUM_NUM FROM khipu_main", con)
    con.close()
    disjoint = get_disjoint_kfg_files(okr_m)
    assert len(disjoint) >= 30
    okr_ids = set(canon_strict(v) for v in okr_m["INVESTIGATOR_NUM"].dropna())
    okr_ids |= set(canon_strict(v) for v in okr_m["orig_inv_num"].dropna())
    for f in disjoint:
        assert canon_strict(f.stem) not in okr_ids, f"Leakage detected: {f.stem}"


def test_phase2a_doseresponse_monotonic_and_significant():
    dr = pd.read_csv(RESULTS_DIR / "phase2a_law5_vertical_doseresponse.csv")
    assert len(dr) == 4
    fracs = dr["frac_top_gt_bot"].tolist()
    # Strictly monotonic climb from 0cm negative control to >=3cm
    assert fracs[0] < fracs[1] < fracs[2] < fracs[3]
    assert fracs[0] < 0.50  # negative control has NO descending preference
    assert fracs[3] > 0.60 and dr.iloc[3]["binom_p"] < 1e-15


def test_phase2a_oos_replication_tables_present():
    l1 = pd.read_csv(RESULTS_DIR / "phase2a_law1_slot_presence_okr_oos.csv")
    oos = l1[l1["cohort"] == "KFG_OOS"]
    assert len(oos) >= 5
    assert int((oos["q_bh"] < 0.05).sum()) >= 4
