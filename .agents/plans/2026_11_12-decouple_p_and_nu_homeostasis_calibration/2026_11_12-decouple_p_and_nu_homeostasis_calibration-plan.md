---
name: "2026_11_12-decouple_p_and_nu_homeostasis_calibration"
description: "Decouple basal division rate (p) from neoplastic division rate (nu) in the ABC-SMC calibration and simulation pipeline to maintain healthy tissue homeostasis"
created_at: "2026-11-12T17:56:30Z"

created_by:
  tool: "Antigravity"
  model:
    name: "Gemini"
    version: "3.8 Flash"
    reasoning_effort: "low"

implemented_by:
  tool: "Antigravity"
  model:
    name: "Gemini"
    version: "3.8 Flash"
    reasoning_effort: "low"

last_implementation_at: "2026-11-12T18:34:00Z"
has_completed_all_phases: false
---

# 🎯 Goal

Decouple the basal cell division rate ($p$) from the neoplastic proliferation rate ($\nu$) in the homeostatic model pipeline, enabling independent calibration where normal tissue division rate balances apoptosis ($p \ge \beta + \mu$) to preserve healthy tissue stability ($N \approx 1000$ cells) while autonomously calibrating neoplastic clonal expansion.

# 👀 Context

- **Mathematical background**:
  - Analysis in `dinamica_texto.txt` proved that normal tissue population obeys $\mathbb{E}[B(t)] = \lambda_B^t B(0)$, with critical boundary $d_h \approx \beta + \mu$.
  - When $p \in [0.047, 0.057]$ and $\nu = 10 \times p$ were strictly tied in `exp08`, the ABC calibration was forced to select $\beta \approx 0.087$ to match clinical penetrance, resulting in $p < \beta$ and causing unintended premature tissue decay.
- **Key files to modify and integrate**:
  - [`experiments/exp08_cs_abc_smc_homeostasis/abc_smc.py`](file:///home/luis/CLionProjects/cellSim/experiments/exp08_cs_abc_smc_homeostasis/abc_smc.py): ABC-SMC calibration engine.
  - [`experiments/exp08_cs_abc_smc_homeostasis/run_homeostasis_features.py`](file:///home/luis/CLionProjects/cellSim/experiments/exp08_cs_abc_smc_homeostasis/run_homeostasis_features.py): Feature evaluation runner for 16,000 virtual tissues.
  - [`experiments/exp08_cs_abc_smc_homeostasis/plot_cell_stage_dynamics.py`](file:///home/luis/CLionProjects/cellSim/experiments/exp08_cs_abc_smc_homeostasis/plot_cell_stage_dynamics.py): Cell stage dynamics visualization.
  - [`experiments/exp08_cs_abc_smc_homeostasis/plot_figure3.py`](file:///home/luis/CLionProjects/cellSim/experiments/exp08_cs_abc_smc_homeostasis/plot_figure3.py): ABC convergence and cumulative risk figure.
  - [`experiments/exp08_cs_abc_smc_homeostasis/plot_figure4_homeostasis.py`](file:///home/luis/CLionProjects/cellSim/experiments/exp08_cs_abc_smc_homeostasis/plot_figure4_homeostasis.py): Publication Figure 4.
  - [`experiments/exp08_cs_abc_smc_homeostasis/plot_sfig4_homeostasis.py`](file:///home/luis/CLionProjects/cellSim/experiments/exp08_cs_abc_smc_homeostasis/plot_sfig4_homeostasis.py): Supplementary Figures S4 and S4b.
- **Architecture and Conventions**:
  - [`AGENTS.md`](file:///home/luis/CLionProjects/cellSim/AGENTS.md): Strict inward dependencies, ADR-0002 agentic cell rules, ADR-0005 Big Bang division.

# 📜 Public Contracts

- **ABC-SMC Parameter Schema (`experiments/exp08_cs_abc_smc_homeostasis/abc_smc.py`)**:
  - Parameter vector $\theta$ expanded from 4 to 5 independent parameters:
    `PARAM_NAMES = ["division_rate", "brca1_rate", "low_delta", "high_delta", "neoplastic_div_rate"]`
  - Priors:
    - `division_rate`: $\mathcal{U}[0.060, 0.140]$ (allows $p \ge \beta + \mu$)
    - `brca1_rate`: $\mathcal{U}[0.001, 0.150]$
    - `low_delta`: $\mathcal{U}[0.060, 0.220]$
    - `high_delta`: $\mathcal{U}[0.120, 0.600]$ (constrained to $\delta_{\rm high} > \delta_{\rm low}$)
    - `neoplastic_div_rate`: $\mathcal{U}[0.150, 0.550]$ (autonomous tumoral proliferation rate)
  - Output CSV format: Preserves `generation_{gen}.csv` columns with all 5 calibrated parameters, weights, distances, and risks at ages 30-80.
- **CLI Commands**:
  - `python3 experiments/exp08_cs_abc_smc_homeostasis/abc_smc.py --particles 100 --generations 5 --n-sim 100 --workers 4`
  - `python3 experiments/exp08_cs_abc_smc_homeostasis/plot_cell_stage_dynamics.py`
  - `python3 experiments/exp08_cs_abc_smc_homeostasis/run_homeostasis_features.py --n-per-particle 80`
  - `python3 experiments/exp08_cs_abc_smc_homeostasis/plot_figure3.py`
  - `python3 experiments/exp08_cs_abc_smc_homeostasis/plot_figure4_homeostasis.py`
  - `python3 experiments/exp08_cs_abc_smc_homeostasis/plot_sfig4_homeostasis.py`

# 🪜 Phases

## Phase 1: Decouple $p$ and $\nu$ in ABC-SMC, recalibrate, and verify healthy tissue preservation

### Description
Modify `abc_smc.py` to decouple `division_rate` ($p$) from `neoplastic_div_rate` ($\nu$), widen the prior for $p$, sample $\nu$ independently, run the 5-parameter ABC-SMC calibration, and verify via cell stage simulation that healthy cell population remains stable ($N \approx 1000$ cells) across early and middle adulthood.

### To-do actions list
- [x] Update `experiments/exp08_cs_abc_smc_homeostasis/abc_smc.py`:
  - Add `neoplastic_div_rate` to `PARAM_NAMES` and `PRIORS`.
  - Expand prior of `division_rate` to $[0.060, 0.140]$ and define `neoplastic_div_rate` prior as $[0.150, 0.550]$.
  - Remove deterministic coupling `nu_div = 10.0 * p_div`.
  - Pass sampled `neoplastic_div_rate` directly to `--neoplastic-div-rate` flag in `run_bootstrapping`.
- [x] Execute ABC-SMC calibration across multiple generations to produce calibrated particles.
- [x] Run `run_stage_trajectories` with the new posterior median parameters ($N=1000$ initial cells, carrying capacity $1000$).
- [x] Execute `plot_cell_stage_dynamics.py` to regenerate `Figure_stage_dynamics.png` and `Figure_stage_dynamics.pdf`, confirming healthy tissue stability ($N \approx 1000$).
- [x] Verify the changes in terms of typechecking, linting and tests using the project's verification command (look it up in the AGENTS.md file or the project configuration). Fix issues if any.
- [x] STOP. Present the changes to the user for review and suggest commit messages (or pull request titles, when the phases are implemented through pull requests). Do NOT proceed to the next phase until the user explicitly asks.

## Phase 2: Simulate emergent clonal features and regenerate manuscript figures

### Description
Run large-scale virtual tissue simulations (`run_homeostasis_features.py`) across the newly decoupled posterior particles, evaluate age-stratified BRCA1 genotype composition and founder lineage complexity, and update the full suite of figures (Figure 3, Figure 4, Figure S4, Figure S4b) reflecting the decoupled homeostatic model.

### To-do actions list
- [ ] Run `experiments/exp08_cs_abc_smc_homeostasis/run_homeostasis_features.py` across calibrated posterior particles to generate updated tissue features CSVs (`Figure4b_genotype_composition_by_age.csv`, `Figure4c_founder_lineages_by_age.csv`, `emergent_16000_tissues.csv`).
- [ ] Update `plot_figure4_homeostasis.py` and `plot_figure4_pairwise_scatters.py` to support the 5-parameter calibrated posterior.
- [ ] Update `plot_sfig4_homeostasis.py` to incorporate $\nu$ in the corner plot and pairwise parameter ratios.
- [ ] Re-run `plot_figure3.py`, `plot_figure4_homeostasis.py`, and `plot_sfig4_homeostasis.py`, exporting updated publication-grade PDFs and PNGs to `results/paper_figures/`.
- [ ] Verify the changes in terms of typechecking, linting and tests using the project's verification command (look it up in the AGENTS.md file or the project configuration). Fix issues if any.
- [ ] STOP. Present the changes to the user for review and suggest commit messages (or pull request titles, when the phases are implemented through pull requests). Do NOT proceed to the next phase until the user explicitly asks.

# ⏭️ Next step

Phase 2 should be completed next.
