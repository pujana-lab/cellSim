---
name: "2026_11_12-abc_pairwise_scatter_plots_figure4"
description: "Generate pairwise parameter scatter plots from ABC-SMC calibration posterior particles to inspect correlations and parameter identifiability"
created_at: "2026-11-12T16:50:00Z"

created_by:
  tool: "Antigravity"
  model:
    name: "Gemini"
    version: "3.8 Flash"
    reasoning_effort: "low"
---

# 🎯 Goal

Generate publication-grade pairwise scatter plots of all calibrated parameters ($p$, $\beta$, $\delta_{\rm low}$, $\delta_{\rm high}$) from the ABC-SMC calibration posterior particles, enabling direct inspection of parameter dependencies, correlations, and potential identifiability through ratios or compensatory relationships.

# 👀 Context

- **Data source**: Posterior particles from Generation 14 of the homeostatic ABC-SMC calibration located at [`experiments/exp08_cs_abc_smc_homeostasis/data/generation_14.csv`](file:///home/luis/CLionProjects/cellSim/experiments/exp08_cs_abc_smc_homeostasis/data/generation_14.csv) (200 particles with columns `division_rate`, `brca1_rate`, `low_delta`, `high_delta`, `weight`).
- **Related scripts**:
  - [`experiments/exp08_cs_abc_smc_homeostasis/plot_figure4_homeostasis.py`](file:///home/luis/CLionProjects/cellSim/experiments/exp08_cs_abc_smc_homeostasis/plot_figure4_homeostasis.py): Generates current Figure 4 with marginal 1D KDEs, CDFs, and correlation heatmap.
  - [`experiments/exp03_cs_abc_smc_calibration/plot_sfig4.py`](file:///home/luis/CLionProjects/cellSim/experiments/exp03_cs_abc_smc_calibration/plot_sfig4.py): Existing corner plot script for the canonical model.
- **Reference documentation**:
  - [`AGENTS.md`](file:///home/luis/CLionProjects/cellSim/AGENTS.md): Architecture, guidelines, and conventions.
  - [`results/paper_figures/README.md`](file:///home/luis/CLionProjects/cellSim/results/paper_figures/README.md): Catalog of manuscript figures.

# 📜 Public Contracts

- **Script CLI**:
  - `experiments/exp08_cs_abc_smc_homeostasis/plot_figure4_pairwise_scatters.py`:
    - Signature / Arguments: `python3 plot_figure4_pairwise_scatters.py [--data-file PATH] [--output-dir PATH]`
- **Generated Visual Artifacts**:
  - `experiments/exp08_cs_abc_smc_homeostasis/Figure4_pairwise_scatters.pdf`: Vector graphics publication format.
  - `experiments/exp08_cs_abc_smc_homeostasis/Figure4_pairwise_scatters.png`: High-resolution (300 DPI) raster image for quick review.
  - `results/paper_figures/Figure4_pairwise_scatters.png` & `.pdf`: Mirrored outputs for the paper figures directory.

# 🪜 Phases

## Phase 1: Pairwise scatter matrix generation and visualization

### Description
Implement a standalone visualization script that loads the calibrated ABC-SMC posterior particles, generates a complete $4 \times 4$ pairwise scatter matrix for $(p, \beta, \delta_{\rm low}, \delta_{\rm high})$ with particle scatter points, 2D density contours (KDE), weighted correlation values ($r$), and marginal histograms/KDEs along the diagonal. Export high-resolution PNG and PDF figures ready to be analyzed.

### To-do actions list
- [x] Create `experiments/exp08_cs_abc_smc_homeostasis/plot_figure4_pairwise_scatters.py` implementing:
  - Particle loading with weights from `generation_14.csv`.
  - Pairwise scatter plots for all 6 parameter combinations:
    - $(p, \beta)$, $(p, \delta_{\rm low})$, $(p, \delta_{\rm high})$, $(\beta, \delta_{\rm low})$, $(\beta, \delta_{\rm high})$, $(\delta_{\rm low}, \delta_{\rm high})$.
  - Weighted bivariate KDE contours and particle scatter points colored/sized by weight.
  - Linear regression trendlines or contour levels showing correlation direction.
  - Annotation of weighted Pearson correlation $r$ in each subplot.
  - Diagonal subplots showing 1D weighted marginal distributions vs uniform prior bounds.
- [x] Execute script to produce `Figure4_pairwise_scatters.png` and `Figure4_pairwise_scatters.pdf` in both `experiments/exp08_cs_abc_smc_homeostasis/` and `results/paper_figures/`.
- [x] Verify that output figures render cleanly and all parameter pairs are clearly readable.
- [x] Verify the changes in terms of typechecking, linting and tests using the project's verification command (look it up in the AGENTS.md file or the project configuration). Fix issues if any.
- [x] STOP. Present the changes to the user for review and suggest commit messages (or pull request titles, when the phases are implemented through pull requests). Do NOT proceed to the next phase until the user explicitly asks.

# ⏭️ Next step

Phase 1 should be completed next.
