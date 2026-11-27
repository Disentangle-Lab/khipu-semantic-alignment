# Pre-Registration: Archive Invariance, Tabular Outliers, and Colonial Schema Alignment in Andean Khipus

> **STATUS:** `FROZEN — HASH RECORDED ON 2026-11-27`  
> **Date frozen:** 2026-11-27  
> **Companion repositories:** `khipu-intra-grammar` (OKR v2.1.0) & `khipu-okr-kfg-concordance` (`n = 611`)

---

## 0. Data Sources, Harmonization & Disjoint Cohorts

- **Primary structural & cord-level data:**
  - `OKR v2.1.0` (`khipu.db`, DOI `10.5281/zenodo.18025748`) for harmonized `CLUSTER_LEVEL = 1` groups, `FULL_COLOR` (`PCORD_FLAG != 1`), `TWIST`, `ATTACHMENT_TYPE`, and intact subsidiary cords (`CORD_LEVEL > 1`).
  - `KFG_DB_Excel/*.xlsx` (from `khipufieldguide/data-science-khipu-article`, **never** `khipu_summary.csv`, to avoid the documented silent drop of 2,313 subsidiary cords) for vertical-position-verified pendant values and cross-repository replication.
  - `master_matched_611.csv` & `deep_dive_fdr_1000.csv` for exact 1-to-1 concordance (`n = 611`) and the pre-established set of **86 doubly-tabular khipus** (`H3a` and `H3b` both BH-FDR `q < 0.05` across the `n = 285` jointly testable OKR khipus).
- **Duplicate & same-object hygiene:** Before any inter-khipu test (Steps 1–3), known duplicate pairs (`UR044`, `AS125/UR247`, `UR066/UR067`, and `khipu-okr-kfg-concordance` `okr_name` collisions) and same-museum-number fragments are flagged; primary analyses are run both **with** and **without** deduplication/same-object collapse.
- **Disjoint Positive Control (Calibration Set):**
  - **Santa Valley (`KH0323–KH0328`, Museo Temple Radicati)** matched against the 1670 colonial tributary census (Medrano & Urton, 2018).
  - Crucially, none of `KH0323–KH0328` belong to the 86 doubly-tabular test set (`KH0323` and `KH0325` pass `H3b` numeric homology at `k = 6`, `p = 0.001`, but have `h3_col_p = 1.0`). Thus Santa Valley is strictly **out-of-sample** with respect to the 86.
- **Test Cohorts:**
  - **Primary cohort:** The `86` doubly-tabular khipus (`H3a ∩ H3b`, BH `q < 0.05`).
  - **Secondary sensitivity cohort (pre-declared, reported separately):** `H3b`-only numeric tabular khipus (`q_val < 0.05`, `q_col >= 0.05`) for numeric-only alignment checks, since known colonial censuses like Santa Valley can be numerically tabular without color-positional homology.

---

## 1. Step 1 — Inter-Khipu Archive Invariance (Controlled for Archive Size & `k`)

### Research Question
Do distinct tabular khipus from the **same archaeological provenance/archive** and with the **same modal column count (`k = h3_m`)** share:
- **(S1a) Column-color sequence** (`FULL_COLOR` modal color per slot `1..k`), and/or
- **(S1b) Column numeric rank profile** (Spearman rank order of column means across slots `1..k`),
beyond what is expected between khipus of the same `k` from *different* archives?

### Pre-Registered Unit & Similarity Metrics
For every khipu `i` with a well-defined `G × k` modal-group matrix (`G >= 2`, `k >= 3`):
- **Color profile** $C_i = (c_{i,1}, \dots, c_{i,k})$: modal `FULL_COLOR` in each column slot `j = 1..k` (ties broken deterministically by first occurrence; also evaluated with full column color-distribution Jensen-Shannon similarity as a robustness check).
- **Numeric rank profile** $R_i = \text{rank}(\bar{x}_{i,1}, \dots, \bar{x}_{i,k})$ where $\bar{x}_{i,j} = \text{mean}(\log(1 + v))$ in column `j`.
- **Pairwise similarity for khipus `(i, i')` with the SAME `k`:**
  - Because khipu reading direction (left-to-right vs. right-to-left) is not universally fixed, compute orientation-invariant similarity:
    $$S_{\text{col}}(i, i') = \max\Big(\frac{1}{k}\sum_{j=1}^{k} \mathbb{1}[c_{i,j} = c_{i',j}],\; \frac{1}{k}\sum_{j=1}^{k} \mathbb{1}[c_{i,j} = c_{i',k+1-j}]\Big)$$
    $$S_{\text{num}}(i, i') = \max\big(\rho_S(R_i, R_{i'}),\; \rho_S(R_i, R_{i'}^{\text{rev}})\big)$$
    (Primary test uses orientation-invariant $\max$; forward-only is reported in the supplement.)

### Controls for Archive Size, `k`, and Cataloger
Preliminary metadata inspection shows strong size skew across `(provenance, k)` cells (e.g., Incahuasi `k=3` has `n=9`; Pachacamac `k=5` has `n=4`, `k=3` has `n=3`; most other provenances have `n=2`). To prevent significance from being an artifact of archive sample size:
1. **Strict `k`-stratified site-label permutation (`N_PERM = 5,000`):** Exclude `Unknown` provenance. Permute provenance labels **only among khipus sharing the exact same `k`**. The test statistic is the **macro-averaged** excess within-site similarity:
   $$T_{\text{macro}} = \frac{1}{|\mathcal{S}|} \sum_{(s,k) \in \mathcal{S}} \Big(\bar{S}_{\text{within}}(s, k) - \bar{S}_{\text{between}}(\cdot, k)\Big)$$
   where $\mathcal{S}$ is the set of valid `(site, k)` cells with `n >= 2`. Because each `(site, k)` cell receives **equal weight `1/|S|`**, large archives (e.g., Incahuasi `k=3` with 36 pairs) cannot drown out small archives (`n=2` with 1 pair).
2. **Matched-size subsampling check:** For every `(site, k)` cell with `n > 2`, repeatedly subsample `n_sub = 2` khipus (1,000 draws) and report the matched-`n` effect size distribution alongside the full-cell estimate.
3. **Leave-one-archive-out & cataloger-stratified control:** Recompute $T_{\text{macro}}$ dropping each top archive (Incahuasi, Pachacamac) one at a time, and run a second permutation restricting label swaps to **same `(k, cataloger)`** strata.

---

## 2. Step 2 — Tabular Taxonomy & Explicit Outlier Accounting (`n = 86`)

### Pre-Registered Table Signature
For each of the 86 doubly-tabular khipus, extract an explicit, interpretable feature vector (combining OKR + rescued KFG Excel values/subsidiaries):
1. **Geometry:** `k = h3_m` (columns), `G = h3_G` (modal rows), coverage fraction `G*k / n_pendants`, aspect ratio `G / k`.
2. **Numeric profile:** column-level monotonicity (Spearman $\rho$ of column means vs. slot index), column zero-rate variance, fraction of groups summing to `×100` (strict & loose), `H5c` cross-group slot-sum flag (`first/last` or `any-row`).
3. **Subsidiary grammar (bug-free):** total subsidiary ratio, and **positional concentration of subsidiaries** across the `k` columns (entropy / chi-square test of whether subsidiaries attach preferentially to specific column slots).
4. **Binary structural axes:** row-level `R/V` purity (are rows purely R or purely V, or mixed?), `S/Z` cord-twist purity, `S/Z` knot-direction balance.

### Family Assignment & Mandatory Outlier Reporting
To avoid forcing every khipu into a cluster:
- **Deterministic structural families** (rule-based, transparent, defined *a priori*):
  - **F1 — Compact parallel registers:** `k ∈ {3, 4}`, `G >= 3`, no `H5c` summary row, `clean_dec_frac >= 0.90`.
  - **F2 — Medium multi-attribute tables:** `k ∈ {5, 6, 7, 8, 9}`, `G >= 3`, no `H5c` summary row, `clean_dec_frac >= 0.90`.
  - **F3 — Self-totaling tables:** Any `k ∈ [3, 12]` with significant `H5c` cross-group slot sums (`p < 0.05`) or strict `×100` group totals (`>= 2` groups).
  - **F4 — Non-clean / hybrid tabular khipus:** `clean_dec_frac < 0.90` (passes H3a+H3b despite non-decimal irregularities).
  - **Outliers (`OUT`):** Any khipu not matching F1–F4 — including extreme arities (`k >= 10`, e.g., `k = 16, 20, 22, 36, 160`), `G = 2` minimal-row tables, or extreme subsidiary/low-coverage structures.
- **Secondary unsupervised check:** HDBSCAN on standardized signature features with **explicit noise (`label == -1`) reporting** (never reassigning noise points to nearest clusters).
- **Mandatory deliverable:** Report exact `n_outliers / 86` (`%`) for both definitions and publish an individual structural dossier for every outlier.

---

## 3. Step 3 — Pre-Registered Colonial Schema Alignment (Table-to-Table)

### Semantic Layer: Structured Colonial Transcriptions (`data/colonial_schemas/`)
Each historical khipu transcription is encoded in JSON with:
- `schema_id`, `source_citation`, `region`, `year`, `domain` (`census`, `mita_labor`, `storehouse_tribute`, `calendar`).
- `k_cols` (number of explicit columns/categories per record) and `m_rows` (number of administrative units/rows, e.g., *ayllus*, *pachakas*, months).
- `value_matrix` (`m_rows × k_cols`, verified against primary scholarly edition, noting any unresolved editorial discrepancies such as Murra 3,938 vs. Julien 4,108).
- `structural_flags`: `has_moiety_split` (*hanan/hurin*), `has_summary_row`, `has_round_100_totals`, `canonical_col_order_desc` (expected relative magnitude rank across columns, derived *strictly from the colonial document*).

Initial pre-registered colonial corpus (target 6–8 verified schemas before running Step 3):
1. **Santa Valley 1670** (Medrano & Urton, 2018) — *Disjoint Positive Control*.
2. **Huánuco / Chupaychu 1549** (Murra, 1982; cross-checked with Julien, 1988).
3. **Visita de Chucuito 1567** (Diez de San Miguel / Murra) — *hanan/hurin* tributary tables.
4. **Visita de Huánuco 1562** (Ortiz de Zúñiga) — village-level household/tribute tables.
5. **Hatun Xauxa 1558–1561** (Murra / Pärssinen) — storehouse & supply accounts.
6. **Sacaca 1579** (Platt / Urton) — tributary/mit'a accounts.
7. **Guaman Poma 1615 / Andean 12-month calendar** — structural reference for `m = 12 / 24`.

### Pre-Registered Alignment Metric (Scale-Invariant, Fixed *a priori*)
To avoid the fragility of exact-value matching across different years/populations while eliminating arbitrary blending weights, a khipu $K$ (`G × k_K`) and a colonial schema $D$ (`m_D × k_D`) are eligible for alignment **only if** $k_K = k_D$ (exact column-arity match; sensitivity check allows $|k_K - k_D| \le 1$ only as a pre-declared secondary analysis).

Let $\bar{p}_K = (\bar{p}_{K,1}, \dots, \bar{p}_{K,k})$ be the khipu's mean column shares (each non-zero row normalized to sum to 1, then averaged across non-zero rows) and $\bar{p}_D$ be the colonial table's mean column shares.

1. **Primary Metric — Orientation-Invariant Jensen-Shannon Similarity ($A_{\text{JSD}} \in [0, 1]$):**
   $$A_{\text{JSD}}(K, D) = 1 - \min\Big(\text{JSD}_2(\bar{p}_K, \bar{p}_D),\; \text{JSD}_2(\bar{p}_K, \bar{p}_D^{\text{rev}})\Big)$$
   where $\text{JSD}_2$ is base-2 Jensen-Shannon divergence (bounded in `[0, 1]`). This single weight-free metric is continuous, scale-invariant, and jointly sensitive to column rank order and relative column shares.
2. **Secondary Robustness Metric (reported separately, never blended):**
   Orientation-invariant Spearman rank correlation (average ranks for ties):
   $$A_{\text{rank}}(K, D) = \max\Big(\rho_S(\bar{p}_K, \bar{p}_D),\; \rho_S(\bar{p}_K, \bar{p}_D^{\text{rev}})\Big)$$
   *(Structural flags such as moiety `R/V` match or `H5c` summary row are NOT folded into $A_{\text{JSD}}$ to avoid researcher degrees of freedom; they are evaluated as an **independent secondary check** on pairs passing FDR, exactly like colour order in `kfg-khipu-analyses` Part B.)*

### Pre-Registered Null Models (`N_PERM = 5,000`, BH-FDR `q < 0.05`)
Because orientation-invariant column permutation yields only $k! / 2$ distinct arrangements (`3` for $k=3$, `12` for $k=4$, `60` for $k=5$, `360` for $k=6$), individual-khipu column-permutation tests are mathematically underpowered to survive BH-FDR for $k \le 5$. We therefore pre-register a resolution-aware dual-null design:

- **Null 2 — Corpus-empirical khipu null (applies individually to ALL $k \ge 3$, `N_PERM = 5,000`):** Compare $A_{\text{JSD}}(K, D)$ against $A_{\text{JSD}}(K', D)$ over random synthetic `G × k` profiles drawn from the empirical pool of khipus sharing $k$ (or within-khipu row-wise column shuffles when fewer than 20 khipus share $k$). Continuous resolution for all $k \ge 3$; controls for generic monotone/Zipf-like column profiles.
- **Null 1 — Column-order permutation null:** Permute the $k$ column labels of $D$ uniformly across all $k!$ permutations (exact enumeration for $k \le 7$; `5,000` random permutations for $k > 7$).
  - For **$k \ge 6$**: evaluated **individually** alongside Null 2 via intersection-union ($p_{\text{both}} = \max(p_{\text{Null1}}, p_{\text{Null2}})$, BH-FDR `q_both < 0.05`).
  - For **$k \in \{3, 4, 5\}$**: individual pairs are evaluated under Null 2 (BH-FDR `q_Null2 < 0.05`), while Null 1 is evaluated at the **aggregate stratum level** (mean $A_{\text{JSD}}$ across all tested khipus of that $k$, or within an archive, against the aggregate permutation distribution, `N_PERM = 5,000`).

### Decision Rule & Stopping Criterion
- **Gate 0 (Positive Control):** Verify behavior on **Santa Valley (`KH0323–KH0328`)** vs. Santa Valley 1670 before running on the 86. Document honestly whether Santa Valley passes or fails; do **not** alter $A_{\text{JSD}}(K, D)$ after inspecting the 86.
- **Confirmatory result:** Only pairs passing the pre-registered FDR criterion (`q_both < 0.05` for $k \ge 6$; `q_Null2 < 0.05` + aggregate Null 1 support for $k \le 5$) are reported as **statistically supported structural matches** (never claimed as literal decipherment without independent archival provenance convergence).
- **Null result commitment:** If `0` pairs among the 86 pass, report that null result directly as evidence that structural table-shape matching alone is insufficient without provenance-matched colonial archives.

---

## 4. Step 4 — Out of Scope (Limitations & Future Work Only)

The **non-decimal / anomalous regime** (`clean_dec_frac < 1.0`, `5,114` non-clean knotted cords in `table3_cord_level_benchmark_611.csv`, and the `42` S/Z spatial calendar candidates) is **explicitly excluded from decoding claims** in this study. Without ground-truth transcriptions for non-decimal/narrative khipus, any semantic attribution remains an unverified hypothesis and will be documented strictly in **Limitations & Future Work**.
