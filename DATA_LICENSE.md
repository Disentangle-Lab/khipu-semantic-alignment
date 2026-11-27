# Data, Figures & Third-Party Redistribution Notice

## Derived Outputs Authored in This Repository

Unless otherwise noted below, all derived CSV tables (`results/`), publication figures (`figures/`), documentation (`README.md`, `RESULTS.md`, `DEVIATIONS.md`, `EXPLORATORY.md`, `PREREGISTRATION.md`), manuscript files (`paper/`), and numerical schema encodings (`data/colonial_schemas/*.json`) authored by Maria Contreras are licensed under the **Creative Commons Attribution 4.0 International License (CC-BY-4.0)**:

https://creativecommons.org/licenses/by/4.0/

Software code (`src/`, `tests/`, `scripts/`) is licensed separately under the **MIT License** (see [`LICENSE`](LICENSE)).

## Third-Party Primary Sources (Not Redistributed in Raw Form)

1. **Open Khipu Repository (OKR v2.1.0):** Not redistributed in this repository (`data/okr/` is gitignored). Download directly from Zenodo (`https://doi.org/10.5281/zenodo.18025748`) under its upstream license.
2. **Khipu Field Guide (KFG) Excel Workbooks:** Not redistributed in this repository (`data/kfg/` is gitignored). Download directly from `https://github.com/khipufieldguide/data-science-khipu-article` under its upstream license.
3. **Textos Andinos Recitations (Pärssinen & Kiviharju, 2004; KFG `textosbook`):** Raw verbatim Spanish transcriptions (`data/textos_andinos_raw/`) are **excluded from this repository and purged from Git history** pending confirmation of redistribution terms with the copyright holders of *Textos Andinos* and the Khipu Field Guide. Users can fetch the source pages directly via `scripts/fetch_data.sh`. Only derived numerical matrices (`data/colonial_schemas/*.json` — aggregate counts, column/row labels, and bibliographic citations without extended paleographic text) are included here for reproducibility; verify redistribution permissions directly with the original publishers before any commercial or bulk reuse.
