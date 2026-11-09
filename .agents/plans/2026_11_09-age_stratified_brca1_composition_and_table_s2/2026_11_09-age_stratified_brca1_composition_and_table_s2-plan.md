---
name: "age_stratified_brca1_composition_and_table_s2"
description: "Replace Figures 4b and 6b with age-stratified BRCA1 composition, replace chi-square test with ordinal trend analysis (Kendall's tau-b), and create Supplementary Table S2 for cell-scale and tissue-scale models."
created_at: "2026-11-09T07:39:52Z"

created_by:
  tool: "Antigravity"
  model:
    name: "Gemini Flash"
    version: "3.8"
    reasoning_effort: "low"

implemented_by:
  tool: "Antigravity"
  model:
    name: "Gemini Flash"
    version: "3.8"
    reasoning_effort: "low"

last_implementation_at: "2026-11-09T08:53:00Z"
has_completed_all_phases: "false"
---

# Plan: Age-Stratified BRCA1 Composition, Ordinal Lineage Trend Analysis, and Supplementary Table S2

## Goal

Enhance the manuscript's clonal and statistical characterization by replacing the single-bar aggregated BRCA1 tumor status figures (Figures 4b and 6b) with age-stratified stacked distributions matching Panels 4c and 6c, replacing the chi-square test on simulated founder counts with ordinal trend association statistics (Kendall's tau-b) and qualitative evaluation, and compiling a comprehensive Supplementary Table S2 documenting all cell-scale and tissue-scale parameters.

## Context

- Developer conventions and architectural constraints are defined in [AGENTS.md](file:///home/luis/CLionProjects/cellSim/AGENTS.md).
- Cell-scale emergent features experiment and figure generation scripts:
  - Data processing: [experiments/exp04_cs_emergent_features/run_emergent_features.py](file:///home/luis/CLionProjects/cellSim/experiments/exp04_cs_emergent_features/run_emergent_features.py)
  - Figure 4 rendering: [experiments/exp04_cs_emergent_features/plot_figure4.py](file:///home/luis/CLionProjects/cellSim/experiments/exp04_cs_emergent_features/plot_figure4.py)
  - Data directory: [experiments/exp04_cs_emergent_features/data/](file:///home/luis/CLionProjects/cellSim/experiments/exp04_cs_emergent_features/data/) containing `emergent_20000_tissues.csv` (20,000 virtual tissues with `onset_year`, `brca1_category`, `tissue_weight`, and `founder_count`).
- Tissue-scale multiobjective cohort experiment and figure generation scripts:
  - Figure 6 rendering: [experiments/exp07_ts_multiobjective_cohort/plot_figure6.py](file:///home/luis/CLionProjects/cellSim/experiments/exp07_ts_multiobjective_cohort/plot_figure6.py)
  - Tissue-scale script: [tissue_scale/scripts/make_figure6_sfig6.py](file:///home/luis/CLionProjects/cellSim/tissue_scale/scripts/make_figure6_sfig6.py)
  - Validation dataset: [experiments/exp07_ts_multiobjective_cohort/data/final_multiobjective_validation_10000.csv](file:///home/luis/CLionProjects/cellSim/experiments/exp07_ts_multiobjective_cohort/data/final_multiobjective_validation_10000.csv) with columns `detect_age`, `brca1_null_fraction`, `major_founders_5pct`.
- Sensitivity and parameter sources across experiments:
  - Global sensitivity: [experiments/exp02_cs_sensitivity_global/run_lhs_sensitivity.py](file:///home/luis/CLionProjects/cellSim/experiments/exp02_cs_sensitivity_global/run_lhs_sensitivity.py)
  - Structural sensitivity: [experiments/exp06_cs_sensitivity_structural/run_structural_sensitivity.py](file:///home/luis/CLionProjects/cellSim/experiments/exp06_cs_sensitivity_structural/run_structural_sensitivity.py)
  - Tissue-scale reference values: [tissue_scale/README.md](file:///home/luis/CLionProjects/cellSim/tissue_scale/README.md)
  - Homogeneous restraint reference: [experiments/exp08_ts_homogeneous_restraint/README.md](file:///home/luis/CLionProjects/cellSim/experiments/exp08_ts_homogeneous_restraint/README.md)
  - Checkpoint relaxation grid: [experiments/exp09_ts_checkpoint_relaxation/run_checkpoint_grid.py](file:///home/luis/CLionProjects/cellSim/experiments/exp09_ts_checkpoint_relaxation/run_checkpoint_grid.py)

## Public Contracts

### Data Transformation & Statistical Contracts
- In `experiments/exp04_cs_emergent_features/run_emergent_features.py`:
  - `summarize_features(raw_df: pd.DataFrame, out_dir: Path)`:
    - Generates `Figure4b_genotype_composition_by_age.csv` with columns: `age_group,exclusively_brca1_pm_pct,mixed_pct,exclusively_brca1_mm_pct`.
    - Updates `Figure4_statistical_tests.csv` to record Kendall's tau-b rank correlation coefficient (`statistic`), two-sided `p_value`, sample size `n_tumors`, and test name `kendall-tau-b`.
- In `experiments/exp07_ts_multiobjective_cohort/plot_figure6.py` and `tissue_scale/scripts/make_figure6_sfig6.py`:
  - `plot_figure_6(...)` and `plot_figure_6_complete(...)`:
    - Bin `brca1_null_fraction` into: `Predominant BRCA1+/-` (<= 0.50), `Mixed` (0.50 to 0.75), `Predominant BRCA1-/-` (> 0.75), stratified by `age_bin` (`<30, 30–39, 40–49, 50–59, 60–69, 70–79, 80`).
    - Compute Kendall's tau-b statistic on `(detect_age, major_founders_5pct)` using `scipy.stats.kendalltau`.

### Visual & Document Artifacts
- Publication Figures:
  - `experiments/exp04_cs_emergent_features/Figure4.{pdf,png,eps}`: Panel 4b shows age-stratified vertical stacked bars matching Panel 4c; Panel 4c reports Kendall's tau-b statistics.
  - `experiments/exp07_ts_multiobjective_cohort/Figure6.{pdf,png,eps}`: Panel 6b shows age-stratified vertical stacked bars matching Panel 6c; Panel 6c reports Kendall's tau-b statistics.
  - Copies synced to `results/paper_figures/Figure4.{pdf,png,eps}` and `results/paper_figures/Figure6.{pdf,png,eps}`.
- Supplementary Tables:
  - `docs/supplementary/Supplementary_Table_S2.md`: Consolidated parameter reference containing Table S2a (Cell-scale) and Table S2b (Tissue-scale) with status categories: `inferred`, `fixed structural assumption`, `exploratory/reference setting`, `derived`, `tissue-scale calibrated`.
  - Machine-readable files: `docs/supplementary/Supplementary_Table_S2a_cell_scale.csv` and `docs/supplementary/Supplementary_Table_S2b_tissue_scale.csv`.

---

## Phases

### Phase 1: Cell-Scale Figure 4 Age Stratification and Statistical Update

#### Description
Update `experiments/exp04_cs_emergent_features/run_emergent_features.py` and `plot_figure4.py` so that Figure 4b displays the proportion of `Exclusively BRCA1+/-`, `Mixed`, and `Exclusively BRCA1-/-` tumors across the 7 age groups (`<30`, `30–39`, `40–49`, `50–59`, `60–69`, `70–79`, `80`), and update Figure 4c to evaluate founder clonal count across age using Kendall's tau-b ordinal trend analysis instead of chi-square.

#### To-do actions list
- [x] In `experiments/exp04_cs_emergent_features/run_emergent_features.py`:
  - Locate `summarize_features(raw_df, out_dir)` (around line 111).
  - Add calculation of weighted genotype percentages per age group using `AGE_ORDER = ["<30", "30-39", "40-49", "50-59", "60-69", "70-79", "80"]`:
    - For each age group `ag`, filter `sub = d[d["age_group"] == ag]`.
    - Compute weighted percentage for each of the three `GENO_ORDER` categories: `Exclusively BRCA1+/-`, `Mixed`, `Exclusively BRCA1-/-`.
    - Save to `out_dir / "Figure4b_genotype_composition_by_age.csv"` with columns: `age_group,exclusively_brca1_pm_pct,mixed_pct,exclusively_brca1_mm_pct`.
  - Replace the `chi2_contingency` block on lines 154-165 with `scipy.stats.kendalltau`:
    - Call `tau_stat, tau_p = stats.kendalltau(d["onset_year"], d["founder_count"])`.
    - Save the test result to `Figure4_statistical_tests.csv` with columns: `test,outcome,statistic,df,p_value,n_tumors` (storing `kendall-tau-b`, `founder count across onset age`, `tau_stat`, `NaN`, `tau_p`, `len(d)`).
  - Run `python3 experiments/exp04_cs_emergent_features/run_emergent_features.py --skip-sim --data-dir experiments/exp04_cs_emergent_features/data --out-dir experiments/exp04_cs_emergent_features/data` to update the CSV files from existing simulation data.
- [x] In `experiments/exp04_cs_emergent_features/plot_figure4.py`:
  - In `render_figure4(data_dir, out_dir)`:
    - Load `df_geno_age = pd.read_csv(data_dir / "Figure4b_genotype_composition_by_age.csv")`.
    - Replace horizontal bar rendering in Panel b (`ax_b` around lines 186-206) with a vertical stacked bar chart across `x = np.arange(len(ages))` using width `0.65`:
      - Bar segments in colors: `#93C5FD` (Exclusively BRCA1+/-), `#FBBF24` (Mixed), `#1E3A8A` (Exclusively BRCA1-/-).
      - Set x-ticks and labels to `["<30", "30–39", "40–49", "50–59", "60–69", "70–79", "80"]` matching Panel c.
      - Set `ax_b.set_ylabel("Tumor-positive tissues (%)", fontsize=9.5)` and `ax_b.set_xlabel("Age at first tumor-positive detection (years)", fontsize=9.5)`.
      - Set `ax_b.set_ylim(0, 100)`.
      - Place legend above or within the panel cleanly.
    - In Panel c (`ax_c` around lines 237-250), update statistical annotation:
      - Read `tau_stat` and `tau_p` from `Figure4_statistical_tests.csv`.
      - Format annotation string as: `rf"Kendall's $\tau_b = {tau_stat:.2f}$ ($P = {tau_p:.1e}$)"` (or `P < 10^{-15}` if infinitesimal).
- [x] Execute `python3 experiments/exp04_cs_emergent_features/plot_figure4.py` to regenerate `Figure4.pdf`, `Figure4.png`, and `Figure4.eps` in `experiments/exp04_cs_emergent_features/`.
- [x] Copy generated Figure 4 files to `results/paper_figures/` (`cp experiments/exp04_cs_emergent_features/Figure4.* results/paper_figures/`).
- [x] Update `experiments/exp04_cs_emergent_features/README.md` to reflect the age-stratified Panel 4b description and replace the chi-square statistic with the Kendall's tau-b rank correlation in the text and summaries.
- [x] Update `results/paper_figures/README.md` to reflect the updated Figure 4 description.
- [x] Verify the changes in terms of typechecking, linting and tests using the project's verification command: `ctest --test-dir build --output-on-failure`. Fix issues if any.
- [ ] STOP. Present the changes to the user for review and suggest commit messages (or pull request titles, when the phases are implemented through pull requests). Do NOT proceed to the next phase until the user explicitly asks.

---

### Phase 2: Tissue-Scale Figure 6 Age Stratification and Statistical Update

#### Description
Update `experiments/exp07_ts_multiobjective_cohort/plot_figure6.py` and `tissue_scale/scripts/make_figure6_sfig6.py` so that Figure 6b displays age-stratified stacked bars for `Predominant BRCA1+/-`, `Mixed`, and `Predominant BRCA1-/-` across the 7 age bins, and update Figure 6c to report Kendall's tau-b rank correlation between `detect_age` and `major_founders_5pct` instead of chi-square.

#### To-do actions list
- [x] In `experiments/exp07_ts_multiobjective_cohort/plot_figure6.py`:
  - Locate Panel b rendering (around line 97).
  - Replace the single horizontal bar calculation with age-stratified cross-tabulation:
    - `pos = final_df[final_df["detect_age"] >= 0].copy()`
    - Define `pos["age_bin"] = pd.cut(pos["detect_age"], bins=[-np.inf, 29, 39, 49, 59, 69, 79, 80], labels=AGE_BIN_LABELS)`
    - Define `pos["brca1_status"] = pd.cut(pos["brca1_null_fraction"], bins=[-np.inf, 0.50, 0.75, np.inf], labels=["Predominant BRCA1+/-", "Mixed", "Predominant BRCA1-/-"])`
    - Compute `tab_b = pd.crosstab(pos["age_bin"], pos["brca1_status"]).reindex(index=AGE_BIN_LABELS, columns=["Predominant BRCA1+/-", "Mixed", "Predominant BRCA1-/-"], fill_value=0)`
    - Compute `pct_b = tab_b.div(tab_b.sum(axis=1), axis=0) * 100.0`
  - Plot vertical stacked bars for Panel b:
    - Colors: `#3498db` (Predominant BRCA1+/-), `#9b59b6` (Mixed), `#e67e22` (Predominant BRCA1-/-).
    - Set x-ticks, x-tick labels (`AGE_BIN_LABELS`), rotation 35 degrees, `ylim(0, 105)`, `ylabel("Tumor-positive tissues (%)")`, `xlabel("Age at first tumor-positive detection (years)")`.
  - In Panel c (lines 124-150):
    - Replace `chi2_contingency` with `tau_stat, tau_p = scipy.stats.kendalltau(pos["detect_age"], pos["major_founders_5pct"])`.
    - Update annotation on `ax_c`: `ax_c.text(0.95, 0.05, rf"Kendall's $\tau_b = {tau_stat:.2f}$ ($P = {tau_p:.1e}$)", transform=ax_c.transAxes, ha="right", fontsize=8.5)`.
- [x] Apply matching changes to `tissue_scale/scripts/make_figure6_sfig6.py` so the standalone tissue-scale package remains completely synchronized.
- [x] Run `python3 experiments/exp07_ts_multiobjective_cohort/plot_figure6.py` and `python3 tissue_scale/scripts/make_figure6_sfig6.py` to regenerate `Figure6.pdf`, `Figure6.png`, and `Figure6.eps`.
- [x] Copy generated Figure 6 files to `results/paper_figures/` (`cp experiments/exp07_ts_multiobjective_cohort/Figure6.* results/paper_figures/`).
- [x] Update `experiments/exp07_ts_multiobjective_cohort/README.md` and `tissue_scale/README.md` to reflect age-stratified Panel 6b and replace the chi-square mention with Kendall's tau-b statistic in Panel 6c.
- [x] Update `results/paper_figures/README.md` to reflect the updated Figure 6 description.
- [x] Verify the changes in terms of typechecking, linting and tests using the project's verification command: `ctest --test-dir build --output-on-failure`. Fix issues if any.
- [ ] STOP. Present the changes to the user for review and suggest commit messages (or pull request titles, when the phases are implemented through pull requests). Do NOT proceed to the next phase until the user explicitly asks.

---

### Phase 3: Compilation and Publication of Supplementary Table S2

#### Description
Compile the unified Supplementary Table S2 partitioned into Table S2a (Cell-scale 500-cell model) and Table S2b (Tissue-scale 10^7-cell model) with exact mathematical symbols, biological meanings, reference values, status classifications (`inferred`, `fixed structural assumption`, `exploratory/reference setting`, `derived`, `tissue-scale calibrated`), biological sources, and sensitivity test references.

#### To-do actions list
- [ ] Ensure `docs/supplementary/` directory exists (`mkdir -p docs/supplementary`).
- [ ] Create `docs/supplementary/Supplementary_Table_S2.md` containing:
  - Header and introductory methodology describing parameter status categories and scale differences.
  - Table S2a (Cell-Scale 500-Cell Agent-Based Model):
    - `brca1_rate` ($\beta$): Inferred 0.0456 yr$^{-1}$ [0.0409, 0.0504]; Status: `inferred`; Basis: ABC-SMC Gen 19; Sensitivity: Exp 02 LHS & PRCC (Fig S1), Exp 01 grid.
    - `low_delta` ($\delta_{\text{low}}$): Inferred 0.1040 [0.0852, 0.1268]; Status: `inferred`; Basis: ABC-SMC; Sensitivity: Exp 02 LHS & RF ~14% (Fig S1).
    - `high_delta` ($\delta_{\text{high}}$): Inferred 0.3875 [0.1397, 0.5853]; Status: `inferred`; Basis: ABC-SMC; Sensitivity: Exp 02 LHS ratio constraint (Fig S1).
    - `neoplastic_div_rate` ($\nu$): Inferred 0.1569 cycle$^{-1}$ [0.1028, 0.2368]; Status: `inferred`; Basis: ABC-SMC; Sensitivity: Exp 02 LHS RF ~15% & locked ablation (Fig S1c), Exp 01.
    - `delta_ratio` ($R_\delta$): 3.65 [1.15, 6.48]; Status: `derived`; Basis: $\delta_{\text{high}} / \delta_{\text{low}}$; Sensitivity: Exp 02 LHS [1.2, 4.0], Fig S4b.
    - `tp53_rate` ($\mu_{\text{TP53}}$): 0.003 yr$^{-1}$; Status: `fixed structural assumption`; Basis: Somatic mutation rate literature & calibration plateau; Sensitivity: Exp 02 LHS RF #1 ~31% (Fig S1), Exp 06 OAT +-10%/+-25% (Fig S2a), Exp 06 bottleneck collapse (Fig S3b).
    - `d1_threshold` ($\theta_{D1}$): 2.0; Status: `fixed structural assumption`; Basis: Two-stage checkpoint barrier (ADR-0002); Sensitivity: Exp 02 LHS (Fig S1), Exp 06 OAT (Fig S2c), Exp 05 interception (Fig 5).
    - `d2_threshold` ($\theta_{D2}$): 5.0; Status: `fixed structural assumption`; Basis: Extrinsic immune evasion threshold; Sensitivity: Exp 02 LHS RF ~9% (Fig S1), Exp 06 OAT (Fig S2b), Exp 05 interception (Fig 5).
    - `division_rate` ($d$): 0.001 cycle$^{-1}$; Status: `fixed structural assumption`; Basis: Resting mammary epithelium turnover; Sensitivity: Exp 02 LHS (Fig S1), Exp 06 OAT (Fig S2d).
    - `kappa_brca1` ($\kappa_{\text{BRCA1}}$): 2.0; Status: `fixed structural assumption`; Basis: Homologous recombination deficiency damage multiplier; Sensitivity: Exp 02 LHS (Fig S1), Exp 06 discrete structural variation (Fig S3a).
    - `age_drift` ($A_{\text{it}}$): 1e-5 yr$^{-1}$; Status: `fixed structural assumption`; Basis: Replicative cellular aging drift; Sensitivity: Exp 02 LHS (Fig S1), Exp 06 OAT (Fig S2f).
    - `tumor_threshold` ($K_{\text{tum}}$): 0.05 (5%); Status: `fixed structural assumption`; Basis: Mammographic detection resolution; Sensitivity: Exp 06 OAT (Fig S2e).
    - `n_cells` ($N_{\text{cells}}$): 500; Status: `fixed structural assumption`; Basis: Terminal ductal lobular unit agent cohort; Sensitivity: Convergence testing across cohorts.
  - Table S2b (Tissue-Scale $10^7$-Cell Population Cohort Model):
    - `N0` ($N_0$): 10,000,000 cells; Status: `fixed structural assumption`; Basis: Intact susceptible human ductal epithelium; Sensitivity: Tissue scaling benchmark.
    - `tumor_positive_cells` ($M_{\text{det}}$): 500,000 cells (5% of $N_0$); Status: `fixed structural assumption`; Basis: 5% clinical tumor mass cutoff; Sensitivity: Preserved invariant across scales.
    - `mu_TP53` ($\mu_{\text{TP53}}$): 2.0076e-5 cycle$^{-1}$; Status: `tissue-scale calibrated`; Basis: Multiobjective calibration in 10,000 tissues; Sensitivity: Exp 08 Pareto frontier (Fig S6), Exp 07 Fig 6.
    - `beta_BRCA1` ($\beta_{\text{BRCA1}}$): 0.00275 yr$^{-1}$; Status: `tissue-scale calibrated`; Basis: Reproducing ~92.5% somatic BRCA1-null tumor composition; Sensitivity: Exp 08 Pareto frontier (Fig S6), Exp 07 Fig 6b.
    - `delta_low` ($\delta_{\text{low}}$): 0.2267; Status: `tissue-scale calibrated`; Basis: Multiobjective calibration; Sensitivity: Exp 08 Fig S6.
    - `delta_high` ($\delta_{\text{high}}$): 1.4928; Status: `tissue-scale calibrated`; Basis: Multiobjective calibration; Sensitivity: Exp 08 Fig S6.
    - `theta_D1` ($\theta_{D1}$): 13.64; Status: `tissue-scale calibrated`; Basis: Cohort accumulation barrier; Sensitivity: Exp 08 Fig S6.
    - `theta_D2_mean` ($\bar{\theta}_{D2}$): 31.79; Status: `tissue-scale calibrated`; Basis: Tissue-restraint barrier mean; Sensitivity: Exp 08 Fig S6.
    - `tau_expansion` ($\tau_T$): 0.65 yr; Status: `tissue-scale calibrated`; Basis: Clinical volumetric doubling times for triple-negative BRCA1 tumors; Sensitivity: Exp 08 Fig S6.
    - `sigma_theta_d2` ($\sigma_{\theta_{D2}}$): 0.38; Status: `tissue-scale calibrated`; Basis: Log-normal inter-individual variation resolving Pareto penetrance trade-off; Sensitivity: Exp 07 Fig 6d, Exp 08 Fig S6.
    - `gamma_driver` ($\gamma$): 0.0 baseline / grid [0, 3e-4]; Status: `exploratory/reference setting`; Basis: Auxiliary driver mutation rate relaxing apoptosis; Sensitivity: Exp 07 Fig 6e, Exp 09 Fig S7.
    - `alpha_bypass` ($p_{\text{bypass}}$): 0.0 baseline / grid [0.25, 1.00]; Status: `exploratory/reference setting`; Basis: Viability survival under BRCA1-first non-canonical order; Sensitivity: Exp 09 Fig S7.
    - Shared baseline rates: $d = 0.001$, $A_{\text{it}} = 10^{-5}$, $\kappa_{\text{BRCA1}} = 2.0$ (Status: `fixed structural assumption`).
  - Cross-Scale Parameter Translation Section explaining the scaling principles between 500-cell and 10^7-cell representations.
- [ ] Export machine-readable CSV files:
  - `docs/supplementary/Supplementary_Table_S2a_cell_scale.csv`
  - `docs/supplementary/Supplementary_Table_S2b_tissue_scale.csv`
- [ ] Reference Supplementary Table S2 in `README.md`, `tissue_scale/README.md`, and `results/paper_figures/README.md`.
- [ ] Verify the changes in terms of typechecking, linting and tests using the project's verification command: `ctest --test-dir build --output-on-failure`. Fix issues if any.
- [ ] STOP. Present the changes to the user for review and suggest commit messages (or pull request titles, when the phases are implemented through pull requests). Do NOT proceed to the next phase until the user explicitly asks.

---

## Next Step

The next step is to implement **Phase 3: Compilation and Publication of Supplementary Table S2**.
