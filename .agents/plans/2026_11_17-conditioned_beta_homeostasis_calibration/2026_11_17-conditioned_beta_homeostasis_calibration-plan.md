---
name: "2026_11_17-conditioned_beta_homeostasis_calibration"
description: "Branch from Thursday baseline d456f6f, condition BRCA1 mutation rate beta on tissue turnover via beta_prior analytical formula to preserve healthy cell dynamics, and execute ABC-SMC calibration"
created_at: "2026-11-17T08:30:30Z"

created_by:
  tool: "Antigravity"
  model:
    name: "Gemini"
    version: "3.8 Flash"
    reasoning_effort: "low"
---

# 🎯 Goal

Branch from Thursday Nov 12 baseline (`d456f6f`), incorporate the clean parameterized execution pipeline, condition the somatic BRCA1 loss rate ($\beta$) on healthy cell turnover ($d$) and TP53 mutation rate ($\mu$) via the analytical formula `beta_prior(mu, d, lo, hi, T=80)`, and calibrate the model using an agile 500-cell setup with 50 particles to reproduce the approved gradual healthy cell decline alongside the Kuchenbaecker penetrance curve and generate the stage dynamics figure.

# 👀 Context

- **Baseline commit**: Commit `d456f6f` ("feat(calibration): decouple basal and neoplastic division rates in homeostatic ABC-SMC") from Thursday Nov 12, which produced the approved 1000-cell stage dynamics figure with gradual baseline cell decline.
- **Mathematical background**:
  - In [`dinamica_texto.txt`](file:///home/luis/CLionProjects/cellSim/dinamica_texto.txt), the expected baseline healthy population evolves as $\mathbb{E}[B(t)] = \lambda_B^t B(0)$, with $\lambda_B = (1 - \beta)(1 - \mu)(1 + d)$.
  - Given $k = (1 - \mu)(1 + d)$, demanding that residual healthy fraction at age $T=80$ satisfies $lo \le \mathbb{E}[B(T)]/B(0) \le hi$ analytically constrains $\beta$ to the interval:
    $$\beta \in \left[1 - \frac{hi^{1/T}}{k},\; 1 - \frac{lo^{1/T}}{k}\right]$$
- **Components to integrate and maintain**:
  - Parameterized runner [`experiments/exp03_cs_abc_smc_calibration/run_pipeline.sh`](file:///home/luis/CLionProjects/cellSim/experiments/exp03_cs_abc_smc_calibration/run_pipeline.sh) with timestamped runs, `--resume-from` support, and JSON metadata.
  - Cell stage dynamics tracker [`cell_scale/app/run_stage_trajectories.cpp`](file:///home/luis/CLionProjects/cellSim/cell_scale/app/run_stage_trajectories.cpp) and visualizer [`experiments/exp03_cs_abc_smc_calibration/plot_cell_stage_dynamics.py`](file:///home/luis/CLionProjects/cellSim/experiments/exp03_cs_abc_smc_calibration/plot_cell_stage_dynamics.py).
  - C++ Domain carrying capacity and decoupled division rates in [`cell_scale/src/domain/tissue/Tissue.cpp`](file:///home/luis/CLionProjects/cellSim/cell_scale/src/domain/tissue/Tissue.cpp) and [`cell_scale/src/domain/cell/AgenticCell.cpp`](file:///home/luis/CLionProjects/cellSim/cell_scale/src/domain/cell/AgenticCell.cpp).
  - Architecture conventions in [`AGENTS.md`](file:///home/luis/CLionProjects/cellSim/AGENTS.md).

# Public Contracts

- **Git Branch**:
  - `experiment/homeostasis-beta-conditioned` branched from `d456f6f`.
- **ABC-SMC Engine (`experiments/exp08_cs_abc_smc_homeostasis/abc_smc.py` or `experiments/exp03_cs_abc_smc_calibration/abc_smc.py`)**:
  - New analytical prior function:
    ```python
    def beta_prior(mu: float, d: float, lo: float, hi: float, T: int = 80) -> Tuple[float, float]:
        k = (1.0 - mu) * (1.0 + d)
        return 1.0 - (hi ** (1.0 / T)) / k, 1.0 - (lo ** (1.0 / T)) / k
    ```
  - CLI flags added:
    - `--healthy-retention-min FLOAT`: Minimum fraction of baseline healthy cells remaining at $T=80$ (default: 0.15).
    - `--healthy-retention-max FLOAT`: Maximum fraction of baseline healthy cells remaining at $T=80$ (default: 0.35).
- **Execution Pipeline (`run_pipeline.sh`)**:
  - Passes `--healthy-retention-min` and `--healthy-retention-max` to `abc_smc.py`.
  - Simulates 1000-cell stage trajectories and exports figures to timestamped run directories.

# 🪜 Phases

## Phase 1: Setup branch, port clean pipeline assets, and implement conditioned `beta_prior`

### Description
Create the target branch `experiment/homeostasis-beta-conditioned` from commit `d456f6f`. Integrate the clean parameterized pipeline assets (`run_pipeline.sh`, `run_stage_trajectories.cpp`, `plot_cell_stage_dynamics.py`, resume capability). Implement the analytical `beta_prior(mu, d, lo, hi)` sampling logic and prior validation in `abc_smc.py`, and run a smoke test to confirm end-to-end functionality.

### To-do actions list
- [x] Create and check out new branch `experiment/homeostasis-beta-conditioned` from `d456f6f`.
- [x] Port the robust parameterized `run_pipeline.sh` script with `--resume-from`, timestamped run directories, and configuration logging.
- [x] Ensure C++ executables `run_bootstrapping` and `run_stage_trajectories` compile cleanly and unit tests pass (`ctest`).
- [x] Implement `beta_prior(mu, d, lo, hi, T=80)` in `abc_smc.py`:
  - When sampling particle priors, first sample `division_rate` ($d$) and `tp53_rate` ($\mu$), then sample `brca1_rate` ($\beta$) uniformly within `beta_prior(mu, d, lo, hi)`.
  - When applying Gaussian perturbation, validate that candidate $\beta$ falls within the dynamically computed `beta_prior(mu_cand, d_cand, lo, hi)`.
- [x] Execute a smoke test (`--smoke-test`) to verify sampling, C++ invocation, and metric calculation without errors.
- [x] Verify the changes in terms of typechecking, linting and tests using the project's verification command (look it up in the AGENTS.md file or the project configuration). Fix issues if any.
- [ ] STOP. Present the changes to the user for review and suggest commit messages (or pull request titles, when the phases are implemented through pull requests). Do NOT proceed to the next phase until the user explicitly asks.

## Phase 2: Execute conditioned ABC-SMC calibration and simulate tissue stage dynamics

### Description
Run the multi-generation ABC-SMC calibration using the conditioned $\beta$ prior with an agile 500-cell setup and a compact particle count ($N=50$ particles, $N_{\rm sim}=200$). Simulate cell stage trajectories for 100 tissues using the resulting MAP posterior parameters, verifying that healthy cells follow the expected gradual decline down to 20-30% at age 80, and render the publication-grade `Figure_stage_dynamics.png`.

### To-do actions list
- [x] Run ABC-SMC calibration with conditioned $\beta$ prior using 600 cells and 50 particles across generations until tolerance converges.
- [x] Extract the MAP parameter vector from the final generation.
- [x] Run `run_stage_trajectories` with MAP parameters for 100 tissues with $N=600$ cells up to age 80.
- [x] Execute `plot_cell_stage_dynamics.py` and generate `Figure_stage_dynamics.png` and `Figure_stage_dynamics.pdf`.
- [x] Verify the changes in terms of typechecking, linting and tests using the project's verification command (look it up in the AGENTS.md file or the project configuration). Fix issues if any.
- [ ] STOP. Present the changes to the user for review and suggest commit messages (or pull request titles, when the phases are implemented through pull requests). Do NOT proceed to the next phase until the user explicitly asks.

## Phase 3: Regenerate publication figures and emergent features suite

### Description
Regenerate Figure 3 (posterior predictive against Kuchenbaecker 2017), Figure S4 (posterior parameter distributions and correlations), and Figure 4 (BRCA1 genotype breakdown and founder lineage complexity) reflecting the conditioned homeostatic model.

### To-do actions list
- [x] Generate Figure 3 (empirical fit and $\varepsilon$ trajectory) and verify that all milestone age penetrances fall within the 95% CI.
- [x] Generate Supplementary Figure S4 and S4b showing the parameter posterior distributions.
- [x] Execute `run_paper_features` or emergent features simulation across posterior particles.
- [x] Render Figure 4 (Panels 4a, 4b, 4c) and synchronize all figures to `results/paper_figures/`.
- [x] Verify the changes in terms of typechecking, linting and tests using the project's verification command (look it up in the AGENTS.md file or the project configuration). Fix issues if any.
- [ ] STOP. Present the changes to the user for review and suggest commit messages (or pull request titles, when the phases are implemented through pull requests). Do NOT proceed to the next phase until the user explicitly asks.

# ⏭️ Next step

Start Phase 1 by creating the new branch `experiment/homeostasis-beta-conditioned` from `d456f6f` and porting the pipeline assets.
