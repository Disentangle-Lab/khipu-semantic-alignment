# Data, Figures & Third-Party Redistribution Notice

## 1. Repository-Authored Derived Outputs (`CC-BY-4.0`)

Unless otherwise noted below, derived CSV tables (`results/`), publication figures (`figures/`), documentation (`README.md`, `RESULTS.md`, `DEVIATIONS.md`, `EXPLORATORY.md`, `PREREGISTRATION.md`), and manuscript files (`paper/`) authored by Maria Contreras in this repository are licensed under the **Creative Commons Attribution 4.0 International License (CC-BY-4.0)**:

https://creativecommons.org/licenses/by/4.0/

Software code authored in this repository (`src/`, `tests/`, `scripts/`) is licensed separately under the **MIT License** (see [`LICENSE`](LICENSE)).

> **Zenodo archive note:** Because Zenodo metadata supports only a single record-level license field (`CC-BY-4.0`), please consult this file and [`LICENSE`](LICENSE) for per-directory terms: code (`src/`, `tests/`, `scripts/`) is MIT-licensed, `data/upstream/` inherits its upstream companion-repository licenses (MIT, documented in Section 2), and `data/colonial_schemas/*.json` are **excluded** from the CC-BY-4.0 grant (Section 3).

---

## 2. Companion Upstream Artifacts Included in `data/upstream/`

The four files in `data/upstream/` are vendored snapshots from companion **Disentangle Lab** repositories so that Phase 1–2 pipelines and unit tests remain self-contained without cloning external repositories. Each file retains the **actual license declared in the root `LICENSE` file of its source repository** (verified directly on GitHub):

| File in `data/upstream/` | Source Repository & Exact Commit | Source Path | Verified Upstream License |
| :--- | :--- | :--- | :--- |
| `controls_h1_h2_h3.py` | [`Disentangle-Lab/khipu-intra-grammar`](https://github.com/Disentangle-Lab/khipu-intra-grammar) @ `eda5008fa8a2c9c9ad72fbc36802cb31ad4fa1f8` (`2026-11-27`) | `src/controls_h1_h2_h3.py` | **MIT License** (`LICENSE`) |
| `deep_dive_fdr_1000.csv` | [`Disentangle-Lab/khipu-intra-grammar`](https://github.com/Disentangle-Lab/khipu-intra-grammar) @ `eda5008fa8a2c9c9ad72fbc36802cb31ad4fa1f8` (`2026-11-27`) | `results/deep_dive_fdr_1000.csv` | **MIT License** (`LICENSE`) |
| `extension_h5_sums.csv` | [`Disentangle-Lab/khipu-intra-grammar`](https://github.com/Disentangle-Lab/khipu-intra-grammar) @ `eda5008fa8a2c9c9ad72fbc36802cb31ad4fa1f8` (`2026-11-27`) | `results/extension_h5_sums.csv` | **MIT License** (`LICENSE`) |
| `master_matched_611.csv` | [`Disentangle-Lab/khipu-okr-kfg-concordance`](https://github.com/Disentangle-Lab/khipu-okr-kfg-concordance) @ `c656a82188cb539a9ba4f442d6628aa0511330ac` (`2026-11-27`) | `results/master_matched_611.csv` | **MIT License** (`LICENSE`) |

---

## 3. Derived Colonial Numerical Schemas (`data/colonial_schemas/*.json`) — Permission Pending, Excluded from CC-BY-4.0

The five JSON files in `data/colonial_schemas/` (`chucuito_1567_moiety_k5.json`, `huaraz_1532_full_k7.json`, `huaraz_1558_demo_k6.json`, `huaraz_1558_full_k7.json`, `sacaca_1551_k22.json`) contain derived numerical matrices (aggregate counts, structural flags, column/row labels, and bibliographic citations without extended paleographic text) extracted from *Textos Andinos* (Pärssinen & Kiviharju, 2004) via KFG `textosbook`.

- **Explicit exclusion:** These JSON files are **NOT** licensed under this repository's `CC-BY-4.0` grant. Formal written redistribution permission from the copyright holders of *Textos Andinos* (Pärssinen & Kiviharju) and the Khipu Field Guide is **pending**.
- Users must verify redistribution permissions directly with the original publishers/authors before any reuse, redistribution, or commercial use beyond scholarly verification of this study.

---

## 4. Third-Party Primary Sources (Not Redistributed in Raw Form)

1. **Open Khipu Repository (OKR v2.1.0):** Not redistributed in this repository (`data/okr/` is gitignored). Download directly from Zenodo (`https://doi.org/10.5281/zenodo.18025748`) via `scripts/fetch_data.sh` (with checksum verification) under its upstream license.
2. **Khipu Field Guide (KFG) Excel Workbooks:** Not redistributed in this repository (`data/kfg/` is gitignored). Download directly from `https://github.com/khipufieldguide/data-science-khipu-article` (pinned to commit `167b8df031dc9b8f1819c7a1148e35836d02016b`) via `scripts/fetch_data.sh` (with SHA-256 verification) under its upstream license.
3. **Raw *Textos Andinos* Recitations (Pärssinen & Kiviharju, 2004; KFG `textosbook`):** Raw verbatim Spanish transcriptions (`data/textos_andinos_raw/`) are **excluded from this repository and purged from Git history** pending written confirmation of redistribution terms. Users can fetch source pages locally via `scripts/fetch_data.sh --with-textos` under upstream terms.
