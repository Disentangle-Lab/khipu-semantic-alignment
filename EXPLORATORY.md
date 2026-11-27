# Exploratory Post-Hoc Extensions (`[EXPLORATORY]`)

> **Date executed:** 2026-11-27 (post-hoc, executed after Phase 1 `PREREGISTRATION.md` was recorded)  
> **Scope:** All analyses below are explicitly designated **`[EXPLORATORY]`** and should be treated as hypothesis-generating or independently replicated candidate findings rather than pre-registered confirmatory tests.

---

## 1. Step 1–2 Hardening Extensions (`src/step1_2_hardening.py`)

1. **Exact intra-row subsidiary slot permutation (`primary_86`, `N_PERM = 5,000`):** Replaces the asymptotic chi-square test mentioned in Step 2 with exact intra-row permutations (`p >= 1/5001 = 0.0002` per specimen). Result: **`27 / 38` (`71.1%`)** on raw counts and **`26 / 34` (`76.5%`)** on binary presence pass BH-FDR `q < 0.05` (exact binomial vs. 5% nominal rate: `p = 2.84e-29`; [`results/step2_subsidiary_perm_5000.csv`](results/step2_subsidiary_perm_5000.csv)).
2. **Outlier functional enrichment (two-sided Fisher exact):** Deterministic geometric outliers (`OUT`, `21 / 86`) are enriched in official calendar khipus (`KFG10 + UR006/KH0242`: `3/21` vs. `1/65`, `OR = 10.67`, `p = 0.0435`) and KFG Part B newly identified paired khipus (`4/21` vs. `2/65`, `OR = 7.41`, `p = 0.0294`), whereas unsupervised HDBSCAN noise (`45 / 86`) shows no enrichment (`OR = 0.29–0.43`, all `p > 0.34`; [`results/step2_outlier_enrichment_fisher.csv`](results/step2_outlier_enrichment_fisher.csv)).
3. **Archive color-complexity & lot audit:** Verifies that within-archive color matches involve complex/dyed codes (Pachacámac `k=6`: `58.3%` complex; `k=9`: `72.2%`; Incahuasi `k=3`: `54.6%` within-site vs. `21.0%` between-site across distinct excavation units), while documenting museum-lot caveats (`KH0469`/`KH0470` in Berlin `VA 63038`, `KH0152`/`KH0156` in consecutive Berlin numbers `VA16635`/`VA16636`; [`results/step1_cells_detailed_audit.csv`](results/step1_cells_detailed_audit.csv)).

---

## 2. Phase 2A — Five Candidate Hierarchical Subsidiary Regularities (`src/phase2a_subsidiary_grammar.py`, `src/p0_audits_and_power.py`)

### Deduplication & Multiple-Testing Control Across the Five Candidate Regularities

To ensure results are not inflated by known or putative duplicate/paired specimens, all five candidate regularities were re-evaluated after removing `9` duplicate/paired specimens (`KH0483`, `KH0304`, `KH0473`, `KH0469`, `KH0470`, `KH0156`, `KH0305`, `KH0503`, `KH0504`; see [`results/phase2a_dedup_effect_sizes.csv`](results/phase2a_dedup_effect_sizes.csv)), and both **Benjamini–Hochberg FDR (`q_bh_5reg`)** and **Holm–Bonferroni (`p_holm_5reg`)** corrections were applied across the five primary regularity tests:

| Candidate Regularity | Deduplicated OKR Sample | Effect Size & 95% CI | Test & Exact / Asymptotic `p` | BH `q` (5 tests) | Holm `p` (5 tests) | Blind KFG-Only OOS (`n=27` khipus, `605` subs) |
| :--- | :--- | :--- | :--- | :---: | :---: | :--- |
| **R1: Slot-locked subsidiary presence** | `39 / 83` khipus (`47.0%`) | `47.0%` (Wilson 95% CI `36.6–57.6%`) passing individual BH `q < 0.05` (`N_PERM = 2,000`, `p_min = 1/2001`) | Exact binomial vs. `5%`: `p = 1.47e-28` | `7.35e-28` | `7.35e-28` | **`5 / 6` (`83.3%`)** BH `q < 0.05` (`p = 1.80e-06`) |
| **R2: 4th-axis subsidiary color** (strict `sub != parent`) | `44 / 65` khipus (`67.7%` pos.) | Median `+3.69 pp` (boot 95% CI `+0.72` to `+7.67`); `r_rb = +0.635`, `Z = +4.35` | Asymptotic Wilcoxon: `p = 1.36e-05` | `1.70e-05` | `2.72e-05` | `3 / 5` pos., median `+18.1 pp` (`QU021 +23.2 pp` non-parent) |
| **R3: 5th-axis subsidiary numeric homology** | `41 / 73` khipus (`56.2%` pos.) | Median excess `R² = +5.13 pp` (boot 95% CI `-2.02` to `+11.67`); `r_rb = +0.428`, `Z = +3.18` | Asymptotic Wilcoxon: `p = 1.47e-03` | `1.47e-03` | `1.47e-03` | **`5 / 5` (`100%`)** pos. (`AS005`: `R²=0.500` vs `0.156`, `p=0.019`) |
| **R4: Hierarchical decay `L1 → L2`** (khipu-level) | `363` khipus (`7,971` pairs) | Median khipu `frac(L2 < L1) = 77.8%` (95% CI `70.0–81.0%`), `r_rb = +0.561`, `Z = +8.93` (pair-level `71.6%`, `r_rb = -0.633`, `Z = -46.6`, asymptotic `p < 1e-300`) | Asymptotic Wilcoxon (khipu): `p = 2.41e-19` | `4.02e-19` | `7.23e-19` | **`72.4%`** (`n=402`, `p = 1.15e-32`; khipu med `87.1%`, `p = 1.94e-04`) |
| **R5: Vertical descending sibling order** (`ΔAP ≥ 3 cm`) | `697 / 1,070` pairs (`137` khipus) | `65.1%` (Wilson 95% CI `62.2–67.9%`) vs. `46.7%` negative control at `ΔAP = 0 cm` (`Δ = +18.5 pp`) | Exact binomial: `p = 2.49e-23` (khipu Wilcoxon `p = 1.06e-03`; 6 catalogers) | `6.23e-23` | `9.97e-23` | KFG omits `ATTACH_POS` (`55.2%` mixed height; consistent with `ΔAP = 0` ctrl) |

*Auxiliary deduplicated subtests (`results/phase2a_dedup_aux_tests.csv`):*
- **R4c (Non-trivial exact part-whole summation, `parent >= 10`, `>= 2` knotted subsidiaries):** `24` exact hits across `n = 1,230` parents vs. `9.92` expected under within-khipu shuffle null (`2.42×` enrichment; **`p = 0.0002`** at `N_PERM = 5,000`, minimum possible `p = 1/5001`).
- **Rejected artifact (cataloger control):** An apparent `R → V` attachment asymmetry (`72.1%` flip) was traced to a single cataloger (`tania`, `1,041 / 1,390` pairs across `6` khipus); all other catalogers showed `> 90%` parent `R/V` conservation.
- **Regional variation across 8 archives (`results/phase2a_regional_subsidiary_dialects.csv`):** Exact mirror rate (`sub == parent`) varies across sites (`χ² = 191.3`, `p = 8.12e-38`: Incahuasi `19.8%` / Huacho `20.4%` vs. Nazca/Ica/Pachacámac `5.6–5.7%`), and khipu-level median attachment height varies significantly (Kruskal–Wallis `H = 21.75`, `p = 0.0028`: Santa Valley `8.5 cm`, Puruchuco `4.8 cm` vs. `0.5–1.5 cm` elsewhere).

---

## 3. Phase 2C — Corpus-Wide `611 × 611` Pair Sweep (`src/phase2c_pair_discovery.py`)

A blind rare-number (`>= 11`, non-multiples of `10`) sweep across OKR v2.1.0 recovered known duplicate/pair controls (`AS125/UR247`, `UR066/UR067`, Puruchuco `UR064/UR068`, Incahuasi `UR255/UR267A` and `UR256/UR267B`) and surfaced two unflagged Berlin museum links (`DUPLICATE_FLAG = 0`) pending archival confirmation:

1. **Putative Berlin Split-Ledger Triptych (3 physical artifacts):**
   - **Piece 1 (Master):** `KH0192` (`AS175`/`UR1175`, Berlin `VA 42518`, Pachacámac, `45` groups of `5` = **Part 1 [Groups 1–21]** + **Part 2 [Groups 22–24, internal Part-3 summary bridge]** + **Part 3 [Groups 25–45]**, median length `65 cm`, `114` subsidiaries).
   - **Piece 2 (Companion matching Part 1):** `KH0470` (`UR233`, Berlin `VA 63038 b`, Huacho, `21` groups of `6` with a `100%`-zero spacer column, median length `40.0 cm`, `18` subsidiaries) matches `UR1175` Part 1 on **`88 / 104` (`84.6%`)** exact cells (`r = 0.9823`, `ρ = 0.9795`, `p <= 1/10001 = 1e-4`; physical length Mann–Whitney `p = 6.21e-28`).
   - **Piece 3 (Companion matching Part 3):** `KH0469` (`UR232`, Berlin `VA 63038(B)`, Huacho, `21` groups of `6` with zero spacer, median length `51.5 cm`, `9` subsidiaries) matches `UR1175` Part 3 on **`73 / 105` (`69.5%`)** exact cells (`r = 0.9446`, `ρ = 0.9565`, `p <= 1/10001 = 1e-4`; physical length Mann–Whitney `p = 1.20e-06`).
   - **Part 2 verification (`results/phase2c_part2_va63038_check.csv`):** `UR1175` Part 2 (Groups 22–24, 15 cords) correlates strongly internally with the column-by-column 7-set sums of Part 3 (`r = 0.990, 0.995, 0.978`, confirming Ascher's Note 5), and **does not appear** in any of the three khipus of Berlin lot `VA 63038` (`UR231/KH0468`, `UR232/KH0469`, `UR233/KH0470`; max `2/15` matches).
2. **Putative unflagged reversed physical duplicate (`KH0121` ↔ `KH0473`):** Reading `KH0473` (`UR237`, Berlin `VA 66834`, provenance unknown) in reverse matches `KH0121` (`AS107`/`UR1107`, Berlin `VA 24370(B)`, Márquez) on **`27 / 29` (`93.1%`)** pendant values (`ρ = 0.9332`) and physical cord lengths (`r = 0.9767`, median `|Δlen| = 1.5 cm`, `66.5` vs. `66.2 cm`). This **suggests** the same Berlin physical artifact was cataloged twice in opposite directions and **suggests** Márquez provenance for `UR237`, pending verification against Ethnologisches Museum Berlin inventory records.
