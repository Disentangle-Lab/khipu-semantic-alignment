# Pre-Registered Archive Invariance, Candidate Subsidiary Regularities, and Putative Split-Ledger Links in Inka Khipus

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24610674.svg)](https://doi.org/10.5281/zenodo.24610674)
[![CI](https://github.com/Disentangle-Lab/khipu-semantic-alignment/actions/workflows/ci.yml/badge.svg)](https://github.com/Disentangle-Lab/khipu-semantic-alignment/actions/workflows/ci.yml)

**Maria Contreras** ([Disentangle Lab](https://github.com/Disentangle-Lab) · [ORCID](https://orcid.org/0009-0003-0749-4971)) · Companion to [khipu-structural-graphs](https://github.com/Disentangle-Lab/khipu-structural-graphs), [khipu-intra-grammar](https://github.com/Disentangle-Lab/khipu-intra-grammar), and [khipu-okr-kfg-concordance](https://github.com/Disentangle-Lab/khipu-okr-kfg-concordance).

This repository evaluates whether surviving Inka tabular khipus share archive-level column schemas (`[PREREG]`), accounts explicitly for structural outliers (`[PREREG]`), tests table-to-table alignment against verbatim highland colonial transcriptions (`[PREREG]`), and reports post-hoc corpus-wide investigations of subsidiary-cord regularities and cross-khipu numerical links (`[EXPLORATORY]`). Full numerical tables and figures live in [`RESULTS.md`](RESULTS.md); every protocol deviation is logged in [`DEVIATIONS.md`](DEVIATIONS.md); post-hoc analyses are segregated in [`EXPLORATORY.md`](EXPLORATORY.md).

---

## Status & Key Takeaways

- **`[PREREG]` Protocol:** [`PREREGISTRATION.md`](PREREGISTRATION.md) · integrity file [`PREREGISTRATION.sha256`](PREREGISTRATION.sha256) (LF, no BOM; same-day freeze/execution caveat documented in [`DEVIATIONS.md`](DEVIATIONS.md)).
- **`[PREREG]` Step 1 (Archive Invariance):** Significant within-archive column invariance controlling for column count `k` and sample size (`S1a` color JSD `p_macro = 0.0034`; `S1b` numeric rank `p_macro = 0.0010`, `N_PERM = 5,000`).
- **`[PREREG]` Step 2 (Outlier Accounting):** `21 / 86` (`24.4%`) fall outside deterministic families (`45 / 86` HDBSCAN noise); individual outlier dossiers are in [`results/step2_outliers_dossier.csv`](results/step2_outliers_dossier.csv).
- **`[PREREG]` Step 3 (Colonial Alignment Null + Power):** `0 / 30` exact-arity (`0 / 70` `|k_K − k_D| == 1`) BH-FDR matches in `primary_86` against five highland *Textos Andinos* schemas, despite `96.7–100%` nominal detection power up to `σ = 0.20` noise ([`results/step3_power_analysis.csv`](results/step3_power_analysis.csv)). **Working hypothesis:** cross-modal alignment is currently bottlenecked by geographic mismatch between predominantly coastal khipus and digitized highland transcriptions rather than absence of tabular syntax.
- **`[EXPLORATORY]` Phase 2A (Subsidiary Regularities):** Correcting KFG `khipu_summary.csv`'s silent drop of `2,313` subsidiary cords, five candidate hierarchical regularities (`R1–R5`) survive specimen deduplication and Holm correction (`p_holm ≤ 1.47e-03`), replicate in direction/magnitude on `36` disjoint KFG-only workbooks, and show regional archive variation across eight sites ([`EXPLORATORY.md`](EXPLORATORY.md)).
- **`[EXPLORATORY]` Phase 2C (Pair Sweep):** A `611 × 611` rare-number sweep recovers known accounting controls and surfaces a putative three-khipu Berlin split-ledger link (`KH0192`/`VA 42518` Parts 1 & 3 ↔ `KH0470`/`VA 63038 b` & `KH0469`/`VA 63038(B)`) plus a putative reversed physical duplicate (`KH0121` ↔ `KH0473`), pending museum verification.
- **Step 4 (Scope Boundary):** Non-decimal and narrative khipus remain strictly out of scope for decoding claims ([`RESULTS.md`](RESULTS.md)).

---

## Quickstart & Data Setup

Raw third-party archives (`OKR v2.1.0`, `KFG Excel workbooks`, and raw *Textos Andinos* texts) are **not redistributed** in Git (see [`DATA_LICENSE.md`](DATA_LICENSE.md)). Fetch them locally before running full pipeline scripts:

```bash
# 1. Install dependencies (Python >= 3.12, uv)
uv sync --group dev --group figures

# 2. Download OKR v2.1.0 (Zenodo) & KFG Excel workbooks into data/okr/ and data/kfg/
bash scripts/fetch_data.sh

# 3. Verify pre-registration checksum (LF, no BOM)
sha256sum -c PREREGISTRATION.sha256

# 4. Run full test suite (or Makefile pipeline)
uv run pytest -q
make all
```

Individual modules run uniformly via `uv run python -m src.<module>`:

```bash
uv run python -m src.step1_archive_invariance
uv run python -m src.step2_taxonomy_outliers
uv run python -m src.step1_2_hardening
uv run python -m src.gate0_santa_valley
uv run python -m src.build_colonial_schemas
uv run python -m src.step3_colonial_alignment
uv run python -m src.make_figures
uv run python -m src.phase2a_subsidiary_grammar
uv run python -m src.phase2c_pair_discovery
uv run python -m src.p0_audits_and_power
```

---

## Test Suite Coverage (`16` tests in `tests/`)

- [`tests/test_metrics.py`](tests/test_metrics.py) (`6` tests): orientation invariance of `S_col` / `S_num`, base-2 JSD bounds `[0, 1]`, macro-averaging invariance to large-cell skew, and deterministic family assignment (`F1–F4` vs. `OUT`).
- [`tests/test_step3_alignment.py`](tests/test_step3_alignment.py) (`4` tests): scale invariance and orientation invariance of `A_JSD`, exact/Monte Carlo Null 1 permutation calibration, and JSON schema matrix validation.
- [`tests/test_phase2a_subsidiary.py`](tests/test_phase2a_subsidiary.py) (`3` tests): strict zero-leakage disjointness between OKR and KFG-only OOS workbooks (`@pytest.mark.requires_raw_data`), monotonic `ΔATTACH_POS` dose-response curve (`0 cm` negative control vs. `≥ 3 cm`), and OOS Law 1 replication table integrity.
- [`tests/test_phase2c_pairs.py`](tests/test_phase2c_pairs.py) (`3` tests): putative Berlin split-ledger cell-match counts and cord-length distinctness, reversed physical duplicate cord-length parity (`KH0121 ↔ KH0473`), and recovery of known positive-control pairs in the `611 × 611` sweep.

---

## License & Citation

- **Code (`src/`, `tests/`, `scripts/`):** MIT License ([`LICENSE`](LICENSE)).
- **Derived tables, figures & manuscript:** CC-BY-4.0 ([`DATA_LICENSE.md`](DATA_LICENSE.md)).
- **Citation:** See [`CITATION.cff`](CITATION.cff) (GitHub *"Cite this repository"* button).
