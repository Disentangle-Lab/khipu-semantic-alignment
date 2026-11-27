"""
data_loader.py — Harmonized loader for OKR v2.1.0 (khipu.db), KFG Excel workbooks,
and upstream concordance / H3 / H5 tables.

Strictly preserves the decoding and schema conventions of:
  - khipu-intra-grammar (controls_h1_h2_h3.py, deep_dive_fdr_h3.py)
  - khipu-okr-kfg-concordance (master_concordance_611.py)
Never uses KFG khipu_summary.csv for subsidiary counts (avoids the documented 2,313-cord drop bug).
"""

from __future__ import annotations

import pathlib
import re
import sqlite3
import warnings
from functools import lru_cache

import numpy as np
import openpyxl
import pandas as pd
from scipy.stats import false_discovery_control

warnings.filterwarnings("ignore")

ROOT = pathlib.Path(__file__).resolve().parent.parent
DATA_DIR = ROOT / "data"
UPSTREAM_DIR = DATA_DIR / "upstream"
KFG_XLSX_DIR = DATA_DIR / "kfg" / "data" / "KFG_DB_Excel"
RESULTS_DIR = ROOT / "results"
FIGURES_DIR = ROOT / "figures"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

GENERIC_MUSEUM_NUMS = {
    "", "1", "sn", "nan", "unnumbered", "no1", "2", "3", "4", "5", "6",
    "4a", "none", "unknown", "na", "nd",
}

KNOWN_DUP_PAIRS = [
    ("AS125", "UR247"),
    ("UR066", "UR067"),
]


def find_okr_db() -> pathlib.Path:
    """Locate khipu.db inside data/okr/ (handles direct file or extracted zip)."""
    okr_dir = DATA_DIR / "okr"
    matches = list(okr_dir.rglob("khipu.db"))
    if not matches:
        matches = list(okr_dir.rglob("*.db")) + list(okr_dir.rglob("*.sqlite"))
    if not matches:
        raise FileNotFoundError(f"Could not find OKR SQLite database under {okr_dir}")
    return matches[0]


def norm_str(val) -> str:
    if val is None or (isinstance(val, float) and np.isnan(val)):
        return ""
    return str(val).strip().strip('"').strip()


def norm_museum_num(val) -> str:
    return re.sub(r"[^a-z0-9]", "", norm_str(val).lower())


def clean_provenance(kfg_prov, okr_prov) -> str:
    for raw in (kfg_prov, okr_prov):
        s = norm_str(raw)
        if s and s.lower() not in {"unknown", "none", "nan", "?", "unprovenienced", "n/a"}:
            return s
    return "Unknown"


def decode_base10_cords(con: sqlite3.Connection) -> pd.DataFrame:
    """Exact base-10 decoder from khipu-intra-grammar/src/controls_h1_h2_h3.py (SQL-accelerated, bit-for-bit identical)."""
    q = """
    SELECT kc.CORD_ID, kc.CLUSTER_ID, kc.ORDINAL AS cl_ord, kc.START_POS,
           COUNT(*) AS n_rows,
           SUM(CASE WHEN UPPER(TRIM(COALESCE(k.TYPE_CODE, ''))) = 'S' THEN 1 ELSE 0 END) AS n_S,
           SUM(CASE WHEN UPPER(TRIM(COALESCE(k.TYPE_CODE, ''))) = 'L' THEN 1 ELSE 0 END) AS n_L,
           SUM(CASE WHEN UPPER(TRIM(COALESCE(k.TYPE_CODE, ''))) = 'E' THEN 1 ELSE 0 END) AS n_E,
           MAX(COALESCE(k.NUM_TURNS, 0)) AS max_turns
    FROM knot_cluster kc
    JOIN knot k ON kc.CLUSTER_ID = k.CLUSTER_ID
    GROUP BY kc.CORD_ID, kc.CLUSTER_ID, kc.ORDINAL, kc.START_POS
    ORDER BY kc.CORD_ID, kc.ORDINAL, kc.START_POS
    """
    rows = con.execute(q).fetchall()
    out = []
    cur_cid = None
    buf = []

    def _flush(cid, clusters):
        n_cl = len(clusters)
        clean = True
        digits = []
        bottom_is_units = False
        for i, (_, _, _, _, n_rows, n_S, n_L, n_E, max_turns) in enumerate(clusters):
            is_bottom = (i == n_cl - 1)
            if n_S == n_rows and 1 <= n_S <= 9:
                digits.append(int(n_S))
            elif n_E == n_rows and n_E == 1 and is_bottom:
                digits.append(1)
                bottom_is_units = True
            elif n_L == n_rows and n_L == 1 and 2 <= max_turns <= 9 and is_bottom:
                digits.append(int(round(max_turns)))
                bottom_is_units = True
            else:
                clean = False
                if n_L > 0 and max_turns > 0:
                    digits.append(int(round(max_turns)))
                    if is_bottom:
                        bottom_is_units = True
                elif n_E > 0:
                    digits.append(int(n_E))
                    if is_bottom:
                        bottom_is_units = True
                else:
                    digits.append(int(n_rows))

        val = 0
        offset = 0 if bottom_is_units else 1
        for pos_from_bottom, d in enumerate(reversed(digits)):
            power = pos_from_bottom + offset
            if power <= 7:
                val += d * (10 ** power)
            else:
                clean = False
        out.append((cid, float(val), int(clean), n_cl))

    for r in rows:
        cid = r[0]
        if cid != cur_cid:
            if cur_cid is not None:
                _flush(cur_cid, buf)
            cur_cid = cid
            buf = [r]
        else:
            buf.append(r)
    if cur_cid is not None:
        _flush(cur_cid, buf)

    return pd.DataFrame(out, columns=["CORD_ID", "cord_val", "is_clean_decimal", "n_knot_clusters"])


def compute_subsidiaries_per_pendant(con: sqlite3.Connection) -> pd.DataFrame:
    """
    Trace every subsidiary cord (CORD_LEVEL > 1) in OKR back to its root level-1 pendant CORD_ID
    via PENDANT_FROM / ATTACHED_TO so each pendant gets an exact bug-free subsidiary count.
    """
    cords_all = pd.read_sql(
        "SELECT CORD_ID, KHIPU_ID, CORD_LEVEL, PENDANT_FROM, ATTACHED_TO FROM cord",
        con,
    )
    lvl1_ids = set(cords_all.loc[cords_all["CORD_LEVEL"] == 1, "CORD_ID"].tolist())
    parent_map = {}
    for row in cords_all.itertuples(index=False):
        if row.CORD_LEVEL and row.CORD_LEVEL > 1:
            # Prefer PENDANT_FROM if it points to another cord; fallback to ATTACHED_TO
            p = row.PENDANT_FROM if (row.PENDANT_FROM in lvl1_ids or row.PENDANT_FROM) else row.ATTACHED_TO
            if p and p != row.CORD_ID:
                parent_map[row.CORD_ID] = p
            elif row.ATTACHED_TO and row.ATTACHED_TO != row.CORD_ID:
                parent_map[row.CORD_ID] = row.ATTACHED_TO

    all_cord_ids = set(cords_all["CORD_ID"].tolist())
    # Map direct parent for multi-level resolution
    direct_parent = {}
    for row in cords_all.itertuples(index=False):
        if row.CORD_LEVEL and row.CORD_LEVEL > 1:
            if row.ATTACHED_TO in all_cord_ids and row.ATTACHED_TO != row.CORD_ID:
                direct_parent[row.CORD_ID] = row.ATTACHED_TO
            elif row.PENDANT_FROM in all_cord_ids and row.PENDANT_FROM != row.CORD_ID:
                direct_parent[row.CORD_ID] = row.PENDANT_FROM

    sub_counts: dict[int, int] = {}
    sub_direct_counts: dict[int, int] = {}
    khipu_sub_totals = (
        cords_all[cords_all["CORD_LEVEL"] > 1].groupby("KHIPU_ID").size().to_dict()
    )

    for cid in direct_parent:
        cur = cid
        seen = {cur}
        steps = 0
        while cur not in lvl1_ids and cur in direct_parent and steps < 12:
            nxt = direct_parent[cur]
            if nxt in seen:
                break
            seen.add(nxt)
            cur = nxt
            steps += 1
        if cur in lvl1_ids:
            sub_counts[cur] = sub_counts.get(cur, 0) + 1
            if steps == 1:
                sub_direct_counts[cur] = sub_direct_counts.get(cur, 0) + 1

    out = pd.DataFrame(
        [
            {"CORD_ID": cid, "n_sub_attached": cnt, "n_sub_direct": sub_direct_counts.get(cid, 0)}
            for cid, cnt in sub_counts.items()
        ]
    )
    return out, khipu_sub_totals


def load_okr_pendants(db_path: pathlib.Path | None = None) -> tuple[pd.DataFrame, dict[int, int]]:
    """Load level-1 pendants with harmonized FULL_COLOR, base-10 values, twist, R/V, knot S/Z, and subsidiaries."""
    if db_path is None:
        db_path = find_okr_db()
    con = sqlite3.connect(str(db_path))
    vals = decode_base10_cords(con)
    sub_df, khipu_sub_totals = compute_subsidiaries_per_pendant(con)

    q = """
    WITH kn AS (
        SELECT CORD_ID,
               COUNT(*) AS n_knots,
               SUM(CASE WHEN UPPER(TRIM(COALESCE(DIRECTION, ''))) = 'S' THEN 1 ELSE 0 END) AS n_knots_S,
               SUM(CASE WHEN UPPER(TRIM(COALESCE(DIRECTION, ''))) = 'Z' THEN 1 ELSE 0 END) AS n_knots_Z
        FROM knot
        GROUP BY CORD_ID
    ),
    col AS (
        SELECT CORD_ID,
               MIN(COALESCE(NULLIF(TRIM(FULL_COLOR), ''), COLOR_CD_1)) AS color_key
        FROM ascher_cord_color
        WHERE COALESCE(PCORD_FLAG, 0) != 1
        GROUP BY CORD_ID
    ),
    col_fallback AS (
        SELECT CORD_ID,
               MIN(COALESCE(NULLIF(TRIM(FULL_COLOR), ''), COLOR_CD_1)) AS color_key_fb
        FROM ascher_cord_color
        GROUP BY CORD_ID
    )
    SELECT
        c.KHIPU_ID, k.OKR_NUM, k.INVESTIGATOR_NUM, k.CREATED_BY, k.MUSEUM_NAME, k.MUSEUM_NUM,
        k.PROVENANCE AS okr_provenance, k.REGION AS okr_region,
        c.CORD_ID, c.CLUSTER_ID, c.CORD_ORDINAL, c.ATTACH_POS,
        UPPER(TRIM(COALESCE(c.TWIST, ''))) AS twist,
        UPPER(TRIM(COALESCE(c.ATTACHMENT_TYPE, ''))) AS att_type,
        UPPER(TRIM(COALESCE(c.TERMINATION, ''))) AS termination,
        COALESCE(col.color_key, col_fallback.color_key_fb) AS color_key,
        COALESCE(kn.n_knots, 0) AS n_knots,
        COALESCE(kn.n_knots_S, 0) AS n_knots_S,
        COALESCE(kn.n_knots_Z, 0) AS n_knots_Z
    FROM cord c
    JOIN khipu_main k ON c.KHIPU_ID = k.KHIPU_ID
    LEFT JOIN col ON c.CORD_ID = col.CORD_ID
    LEFT JOIN col_fallback ON c.CORD_ID = col_fallback.CORD_ID
    LEFT JOIN kn ON c.CORD_ID = kn.CORD_ID
    WHERE c.CORD_LEVEL = 1
    ORDER BY c.KHIPU_ID, c.CORD_ORDINAL, c.ATTACH_POS, c.CORD_ID
    """
    df = pd.read_sql(q, con)
    con.close()

    df = df.merge(vals, on="CORD_ID", how="left")
    if len(sub_df) > 0:
        df = df.merge(sub_df, on="CORD_ID", how="left")
    else:
        df["n_sub_attached"] = 0
        df["n_sub_direct"] = 0

    df["cord_val"] = df["cord_val"].fillna(0.0)
    df["is_clean_decimal"] = df["is_clean_decimal"].fillna(1).astype(int)
    df["n_sub_attached"] = df["n_sub_attached"].fillna(0).astype(int)
    df["n_sub_direct"] = df["n_sub_direct"].fillna(0).astype(int)
    df["cataloger"] = (
        df["INVESTIGATOR_NUM"].fillna("").astype(str).str.extract(r"^([A-Za-z]{2})")[0].str.upper().fillna("NA")
    )
    # Normalize color strings (strip spaces, uppercase)
    df["color_norm"] = (
        df["color_key"].fillna("").astype(str).str.strip().str.upper().str.replace(r"\s+", "", regex=True)
    )
    df.loc[df["color_norm"] == "", "color_norm"] = "__MISSING__"
    return df, khipu_sub_totals


@lru_cache(maxsize=1024)
def load_kfg_excel_pendants(legal_filename: str) -> pd.DataFrame | None:
    """Load KFG Excel workbook level-1 pendants with KFG vertical-position values, clusters, and terminations."""
    if not legal_filename:
        return None
    path = KFG_XLSX_DIR / f"{legal_filename}.xlsx"
    if not path.exists():
        return None
    wb = openpyxl.load_workbook(path, read_only=True, data_only=True)
    clusters = []
    if "Clusters" in wb.sheetnames:
        for (r,) in wb["Clusters"].iter_rows(values_only=True, max_col=1):
            if not r or str(r).startswith("!"):
                continue
            m = re.search(r"\((\d+)(?:-(\d+))?\)", str(r))
            if m:
                clusters.append((int(m.group(1)), int(m.group(2) or m.group(1))))

    if "Cords" not in wb.sheetnames:
        return None
    it = wb["Cords"].iter_rows(values_only=True)
    hdr = None
    for r in it:
        if r and r[0] == "Cord_Name":
            hdr = {str(x).strip(): i for i, x in enumerate(r) if x is not None}
            break
    if hdr is None:
        return None

    rows = []
    for r in it:
        if not r or r[0] is None:
            continue
        cn = str(r[0]).strip()
        if cn.startswith("!"):
            continue
        if cn.count("s") > 0:
            continue  # pendant only
        m_idx = re.search(r"p(\d+)", cn)
        if not m_idx:
            continue
        idx = int(m_idx.group(1))
        raw_v = r[hdr["Value"]] if "Value" in hdr and hdr["Value"] < len(r) else None
        try:
            v = int(float(raw_v)) if raw_v not in (None, "") else 0
        except (TypeError, ValueError):
            v = 0
        col = str(r[hdr["Color"]] or "").strip().upper() if "Color" in hdr and hdr["Color"] < len(r) else ""
        term = str(r[hdr["Termination"]] or "").strip().upper() if "Termination" in hdr and hdr["Termination"] < len(r) else ""
        rows.append((idx, cn, v, col, term))

    if not rows:
        return None
    p = pd.DataFrame(rows, columns=["p_idx", "cord_name", "kfg_val", "kfg_color", "kfg_term"])
    p = p.sort_values("p_idx").reset_index(drop=True)
    p["kfg_cluster"] = -1
    for ci, (a, b) in enumerate(clusters):
        p.loc[p["p_idx"].between(a, b), "kfg_cluster"] = ci
    return p


def load_cohorts() -> pd.DataFrame:
    """
    Build the specimen-level cohort table uniting:
      - deep_dive_fdr_1000.csv (512 OKR khipus, H3a/H3b FDR -> exact 86 doubly-tabular)
      - master_matched_611.csv (KFG concordance, provenance, legal_filename, museum_num)
      - extension_h5_sums.csv (H5b/H5c exact sum stats)
    Adds deduplication / same-object collapse flags.
    """
    dd = pd.read_csv(UPSTREAM_DIR / "deep_dive_fdr_1000.csv")
    mm = pd.read_csv(UPSTREAM_DIR / "master_matched_611.csv")
    h5 = pd.read_csv(UPSTREAM_DIR / "extension_h5_sums.csv")

    # Select clean non-colliding columns from master_matched_611
    mm_sub = mm[
        [
            "KHIPU_ID", "INVESTIGATOR_NUM", "MUSEUM_NUM", "REGION", "PROVENANCE",
            "kfg_name", "okr_name", "legal_filename", "nickname",
            "provenance", "region", "okr_num_sub", "okr_num_pendant",
        ]
    ].rename(
        columns={
            "PROVENANCE": "okr_prov",
            "REGION": "okr_reg",
            "provenance": "kfg_prov",
            "region": "kfg_reg",
        }
    )

    df = dd.merge(mm_sub, on="KHIPU_ID", how="left")
    df = df.merge(
        h5[["KHIPU_ID", "h5b_hits", "h5b_p", "h5c_G", "h5c_m", "h5c_hits", "h5c_p", "h5c_any_hits", "h5c_any_p"]],
        on="KHIPU_ID",
        how="left",
    )

    # Exact reproduction of H3a (n=299), H3b (n=326), and joint doubly-tabular (n=285 -> 86)
    df["q_col_all"] = np.nan
    m_col = df["h3_col_p"].notna()
    df.loc[m_col, "q_col_all"] = false_discovery_control(df.loc[m_col, "h3_col_p"].values, method="bh")
    df["h3a_sig"] = df["q_col_all"] < 0.05

    df["q_val_all"] = np.nan
    m_val = df["h3_val_p"].notna()
    df.loc[m_val, "q_val_all"] = false_discovery_control(df.loc[m_val, "h3_val_p"].values, method="bh")
    df["h3b_sig"] = df["q_val_all"] < 0.05

    df["q_col_both"] = np.nan
    df["q_val_both"] = np.nan
    m_both = m_col & m_val
    df.loc[m_both, "q_col_both"] = false_discovery_control(df.loc[m_both, "h3_col_p"].values, method="bh")
    df.loc[m_both, "q_val_both"] = false_discovery_control(df.loc[m_both, "h3_val_p"].values, method="bh")
    df["doubly_tabular"] = (df["q_col_both"] < 0.05) & (df["q_val_both"] < 0.05)
    df["h3b_only"] = df["h3b_sig"] & (~df["doubly_tabular"])

    df["provenance_clean"] = [
        clean_provenance(kp, op) for kp, op in zip(df["kfg_prov"], df["okr_prov"])
    ]
    df["museum_clean"] = df["museum"].map(norm_str)
    df["museum_num_norm"] = df["MUSEUM_NUM"].map(norm_museum_num)

    # Deduplication & same-object flag (keep specimen with most pendants per duplicate/museum-object)
    drop_ids = set()
    kfg_to_row = {str(r.kfg_name).strip(): r for r in df.itertuples(index=False) if pd.notna(r.kfg_name)}
    for a, b in KNOWN_DUP_PAIRS:
        if a in kfg_to_row and b in kfg_to_row:
            ra, rb = kfg_to_row[a], kfg_to_row[b]
            loser = rb.KHIPU_ID if ra.n_pendants >= rb.n_pendants else ra.KHIPU_ID
            drop_ids.add(loser)
    if "UR044" in kfg_to_row:
        drop_ids.add(kfg_to_row["UR044"].KHIPU_ID)

    # Same museum + same non-generic museum_num -> same physical object / sibling bundle
    same_obj_drop = set()
    valid_obj = df[(~df["museum_num_norm"].isin(GENERIC_MUSEUM_NUMS)) & (df["museum_num_norm"].str.len() >= 3)]
    for (_, _), grp in valid_obj.groupby(["museum_clean", "museum_num_norm"]):
        if len(grp) > 1:
            ordered = grp.sort_values(["n_pendants", "KHIPU_ID"], ascending=[False, True])
            same_obj_drop.update(ordered["KHIPU_ID"].iloc[1:].tolist())

    df["is_known_dup_drop"] = df["KHIPU_ID"].isin(drop_ids)
    df["is_same_obj_drop"] = df["KHIPU_ID"].isin(same_obj_drop)
    df["keep_dedup"] = (~df["is_known_dup_drop"]) & (~df["is_same_obj_drop"])
    return df


def extract_modal_table(g: pd.DataFrame) -> dict | None:
    """
    Extract the exact G x k modal equal-sized group matrix for a khipu's level-1 pendants,
    matching deep_dive_fdr_h3.py (valid_sizes >= 3, modal_m = value_counts().index[0], G >= 2).
    """
    cl_sizes = g.groupby("CLUSTER_ID", sort=False).size()
    valid_sizes = cl_sizes[cl_sizes >= 3]
    if len(valid_sizes) < 2:
        return None
    sc = valid_sizes.value_counts()
    modal_m = int(sc.index[0])
    if sc.iloc[0] < 2:
        return None
    hom_cids = valid_sizes[valid_sizes == modal_m].index.tolist()

    mat_col = []
    mat_val = []
    mat_att = []
    mat_twist = []
    mat_sub = []
    mat_ks = []
    mat_kz = []
    for cid in hom_cids:
        sub = g[g["CLUSTER_ID"] == cid]
        mat_col.append(sub["color_norm"].values)
        mat_val.append(sub["cord_val"].values.astype(float))
        mat_att.append(sub["att_type"].values)
        mat_twist.append(sub["twist"].values)
        mat_sub.append(sub["n_sub_attached"].values.astype(int))
        mat_ks.append(sub["n_knots_S"].values.astype(int))
        mat_kz.append(sub["n_knots_Z"].values.astype(int))

    return {
        "G": len(hom_cids),
        "k": modal_m,
        "hom_cids": hom_cids,
        "mat_col": np.array(mat_col, dtype=object),
        "mat_val": np.array(mat_val, dtype=float),
        "mat_att": np.array(mat_att, dtype=object),
        "mat_twist": np.array(mat_twist, dtype=object),
        "mat_sub": np.array(mat_sub, dtype=int),
        "mat_ks": np.array(mat_ks, dtype=int),
        "mat_kz": np.array(mat_kz, dtype=int),
    }
