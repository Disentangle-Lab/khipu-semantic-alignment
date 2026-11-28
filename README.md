# Pre-Specified Archive Invariance, Candidate Subsidiary Regularities, and Putative Split-Ledger Links in Inka Khipus

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24610674.svg)](https://doi.org/10.5281/zenodo.24610674)
[![CI](https://github.com/Disentangle-Lab/khipu-semantic-alignment/actions/workflows/ci.yml/badge.svg)](https://github.com/Disentangle-Lab/khipu-semantic-alignment/actions/workflows/ci.yml)

**Maria Contreras** ([Disentangle Lab](https://github.com/Disentangle-Lab) · [ORCID](https://orcid.org/0009-0003-0749-4971)) · Companion to [khipu-structural-graphs](https://github.com/Disentangle-Lab/khipu-structural-graphs), [khipu-intra-grammar](https://github.com/Disentangle-Lab/khipu-intra-grammar), and [khipu-okr-kfg-concordance](https://github.com/Disentangle-Lab/khipu-okr-kfg-concordance).

This repository evaluates whether surviving Inka tabular khipus share archive-level column schemas (`[PREREG]`, internally pre-specified), accounts explicitly for structural outliers (`[PREREG]`), tests table-to-table alignment against verbatim highland colonial transcriptions (`[PREREG]`), and reports post-hoc corpus-wide investigations of subsidiary-cord regularities and cross-khipu numerical links (`[EXPLORATORY]`). Full numerical tables and figures live in [`RESULTS.md`](RESULTS.md); every protocol deviation and filename rename is logged in [`DEVIATIONS.md`](DEVIATIONS.md); post-hoc analyses are segregated in [`EXPLORATORY.md`](EXPLORATORY.md); the preprint manuscript draft is in [`paper/main.tex`](paper/main.tex).

---

## Status & Key Takeaways

- **`[PREREG]` Protocol (Internally Pre-Specified):** [`PREREGISTRATION.md`](PREREGISTRATION.md) · integrity file [`PREREGISTRATION.sha256`](PREREGISTRATION.sha256) (LF, no BOM; same-day freeze/execution caveat documented in [`DEVIATIONS.md`](DEVIATIONS.md)).
- **`[PREREG]` Step 1 (Archive Invariance):** Controlling for column count `k` and sample size (`N_PERM = 5,000`), **primary** `S1a` modal column-color match yields `T_macro = +0.0934` (`p_macro = 0.0128`; drops to `p = 0.2853` when Pachacámac is excluded), while the full-distribution color JSD robustness metric yields `T_macro = +0.1265` (`p_macro = 0.0034`; `p = 0.0326` without Pachacámac). **Primary** `S1b` numeric column rank yields `T_macro = +0.2355` (`p_macro = 0.0010`; `p = 0.0004` without Pachacámac).
- **`[PREREG]` Step 2 (Outlier Accounting):** `21 / 86` (`24.4%`) fall outside deterministic families (`45 / 86` HDBSCAN noise); individual outlier dossiers are in [`results/step2_outliers_dossier.csv`](results/step2_outliers_dossier.csv).
- **`[PREREG]` Step 3 (Colonial Alignment Null + Power):** `0 / 30` exact-arity (`0 / 70` `|k_K − k_D| == 1`) BH-FDR matches in `primary_86` against five highland *Textos Andinos* schemas (3 effective archival traditions: Chucuito, Huaraz, Sacaca). At `20%` multiplicative noise (`σ = 0.20`), **BH-FDR (`m = 30`) detection power** is **`70.8%` (`k = 5`)**, **`59.2%` (`k = 7`)**, and **`0%` for an isolated `k = 6` hit** (due to the `2/720 = 0.00278` Null 1 discreteness floor), alongside `96.7–100%` nominal `p < 0.05` power ([`results/step3_power_analysis.csv`](results/step3_power_analysis.csv)). **Working hypothesis:** cross-modal alignment is currently bottlenecked by geographic mismatch between predominantly coastal khipus and digitized highland transcriptions rather than absence of tabular syntax.
- **`[EXPLORATORY]` Phase 2A (Candidate Subsidiary Regularities):** Correcting KFG `khipu_summary.csv`'s silent drop of `2,313` subsidiary cords, five candidate hierarchical regularities (`R1–R5`) survive specimen deduplication and Holm correction across the five reported khipu-level tests (`p_holm ≤ 2.12e-03`, including khipu-level `R5` Wilcoxon `p = 1.06e-03` across `137` khipus; note that post-hoc Holm over 5 selected tests underestimates broader exploratory multiplicity and `R3`'s median bootstrap 95% CI `-2.02` to `+11.67 pp` crosses 0), replicate in direction/magnitude on `36` disjoint KFG-only workbooks, and show regional archive variation across eight sites ([`EXPLORATORY.md`](EXPLORATORY.md)).
- **`[EXPLORATORY]` Phase 2C (Descriptive Pair Sweep):** A corpus-wide `611 × 611` rare-number sweep (`~186k` pairs; exact pre-filter `n_rare_dist ≥ 5` and [`rare_dice ≥ 0.25` or `n_rare_dist ≥ 10`] yielding `538` candidate pairs in [`results/phase2c_corpus_pair_sweep.csv`](results/phase2c_corpus_pair_sweep.csv)) recovers known controls and ranks two putative Berlin links at **`#1 / 538`** (`KH0192`/`VA 42518` Part 1 ↔ `KH0470`/`VA 63038 b`, with Part 3 ↔ `KH0469`/`VA 63038(B)` identified via targeted lot inspection) and **`#3 / 538`** (putative reversed physical duplicate `KH0121`/`VA 24370(B)` ↔ `KH0473`/`VA 66834`), pending museum verification. Reported `p ≤ 1e-4` values are conditional post-selection permutations uncorrected for sweep-wide multiplicity; a full-sweep empirical null is planned for the companion OKR↔KFG / `v1.1.0` study. **Interactive cell viewer:** [`docs/berlin_triptych_viewer.html`](docs/berlin_triptych_viewer.html) (live GitHub Pages view: [Berlin Triptych Viewer](https://disentangle-lab.github.io/khipu-semantic-alignment/berlin_triptych_viewer.html); rebuilt reproducibly via `uv run python scripts/build_berlin_viewer.py`).
- **Step 4 (Scope Boundary):** Non-decimal and narrative khipus remain strictly out of scope for decoding claims ([`RESULTS.md`](RESULTS.md)).

---

## Quickstart & Data Setup

Raw third-party archives (`OKR v2.1.0`, pinned `KFG` commit, and raw *Textos Andinos* texts) are **not redistributed** in Git (see [`DATA_LICENSE.md`](DATA_LICENSE.md)). Fetch and verify them locally before running full pipeline scripts:

```bash
# 1. Install dependencies (Python >= 3.12, uv)
uv sync --group dev --group figures

# 2. Download OKR v2.1.0 (Zenodo) & pinned KFG commit with SHA-256 verification
bash scripts/fetch_data.sh

# 3. Verify pre-specification checksum (LF, no BOM)
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
- [`tests/test_phase2a_subsidiary.py`](tests/test_phase2a_subsidiary.py) (`3` tests): strict zero-leakage disjointness between OKR and KFG-only OOS workbooks (`@pytest.mark.requires_raw_data`), monotonic `ΔATTACH_POS` dose-response curve (`0 cm` negative control vs. `≥ 3 cm`), and OOS `R1` replication table integrity.
- [`tests/test_phase2c_pairs.py`](tests/test_phase2c_pairs.py) (`3` tests): putative Berlin split-ledger cell-match counts and cord-length distinctness, putative reversed physical duplicate cord-length parity (`KH0121 ↔ KH0473`), and recovery of known positive-control pairs in the `611 × 611` sweep.

---

## License & Citation

- **Code (`src/`, `tests/`, `scripts/`):** MIT License ([`LICENSE`](LICENSE)).
- **Derived tables, figures & manuscript:** CC-BY-4.0 ([`DATA_LICENSE.md`](DATA_LICENSE.md); see explicit separate notices for `data/upstream/` and `data/colonial_schemas/*.json`).
- **Citation:** See [`CITATION.cff`](CITATION.cff) (GitHub *"Cite this repository"* button).
