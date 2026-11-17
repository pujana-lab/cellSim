---
name: "2026_11_11-abc_smc_p0052_homeostasis_calibration"
description: "Implement healthy tissue homeostasis with a 1000-cell carrying capacity cap, tied tumoral division rate (nu = 10 * p), broadened mutation priors, and a calibrated ABC-SMC pipeline against Kuchenbaecker (2017)"
created_at: "2026-11-11T09:46:00Z"

created_by:
  tool: "Antigravity CLI"
  model:
    name: "Gemini"
    version: "3.8 Flash"
    reasoning_effort: "low"
---

# ABC-SMC Calibration with Healthy Tissue Homeostasis (1000 Cells & Division Coupling)

## Goal
Implement a carrying-capacity homeostasis mechanism for normal/healthy cells capping the population at 1000 cells while allowing autonomous neoplastic cell expansion. Couple the tumoral division rate as strictly 10x the basal division rate ($\nu = 10 \times p$, with $p \in [0.047, 0.057]$), broaden somatic mutation priors, and configure an ABC-SMC calibration pipeline against Kuchenbaecker (2017).

## Context
Recent mathematical analysis of tissue dynamics (`brca1_onset_dinamica_tejido.pdf`) demonstrated that the previous 500-cell parametrization ($d = 0.001$, $\beta = 0.0457$) suffered from critical tissue decay due to an inversion of scales ($\beta / d \approx 46$), where normal tissue halved every 14 years. To establish biological realism:
1. Normal cells must maintain homeostatic equilibrium ($N_{\text{normal}} \approx 1000$).
2. Divisions must be more frequent than mutations ($d \gg \beta$).
3. Tumoral cells must divide $\approx 10\times$ faster than normal cells under Big Bang mode.
4. Normal cells must undergo contact inhibition when the normal niche reaches 1000 cells, whereas tumoral cells must remain unrestricted.

Key files:
- [ICell.h](file:///home/luis/CLionProjects/cellSim/cell_scale/src/domain/cell/ICell.h)
- [AgenticCell.h](file:///home/luis/CLionProjects/cellSim/cell_scale/src/domain/cell/AgenticCell.h) / [AgenticCell.cpp](file:///home/luis/CLionProjects/cellSim/cell_scale/src/domain/cell/AgenticCell.cpp)
- [Tissue.h](file:///home/luis/CLionProjects/cellSim/cell_scale/src/domain/tissue/Tissue.h) / [Tissue.cpp](file:///home/luis/CLionProjects/cellSim/cell_scale/src/domain/tissue/Tissue.cpp)
- [run_bootstrapping.cpp](file:///home/luis/CLionProjects/cellSim/cell_scale/app/run_bootstrapping.cpp)
- [abc_smc.py](file:///home/luis/CLionProjects/cellSim/experiments/exp03_cs_abc_smc_calibration/abc_smc.py)

---

## Phases

### Phase 1: Git Branch Creation and C++ Domain Homeostasis

#### Description
Create the isolated git branch for this experiment. Modify the cell and tissue domain model to implement contact inhibition on healthy cells when the healthy population reaches 1000, while preserving unconditional proliferation for neoplastic/tumoral cells. Expose the carrying capacity flag in the CLI bootstrapping runner and write unit tests for the homeostatic equilibrium.

#### Public Contracts
- **`ICell.h`**
  - Add `virtual void setDivisionInhibited(bool inhibited) {}`
  - Add `[[nodiscard]] virtual bool isDivisionInhibited() const { return false; }`
- **`AgenticCell.h` & `AgenticCell.cpp`**
  - Add private member `bool division_inhibited_ = false;`
  - Override `setDivisionInhibited(bool inhibited)` and `isDivisionInhibited()`
  - In `phase5_Cytokinesis()`: Only call `attemptDivision()` if `is_neoplastic_ || !division_inhibited_`.
- **`Tissue.h` & `Tissue.cpp`**
  - Add private member `std::size_t max_healthy_cells_ = 1000;`
  - Add public methods:
    - `void setMaxHealthyCells(std::size_t cap);`
    - `[[nodiscard]] std::size_t getMaxHealthyCells() const;`
  - In `Tissue::live()`:
    - Count living non-neoplastic cells (`!cell->isNeoplastic()`).
    - Set `cell->setDivisionInhibited(living_healthy_count >= max_healthy_cells_)` for non-neoplastic cells.
    - When appending `new_daughters`: Discard non-tumoral daughter if current healthy count reaches `max_healthy_cells_`; always accept tumoral daughters.
- **`run_bootstrapping.cpp`**
  - Add `--max-healthy-cells <int>` flag (default 1000).
  - Forward `max_healthy_cells` to `Tissue`.
- **Unit Test**: `cell_scale/tests/TissueHomeostasisTest.cpp`
  - Test healthy cells do not exceed 1000 even with high division rate $p = 0.052$.
  - Test tumoral cells can divide and expand beyond 1000.

#### To-do Actions List
- [ ] Create and checkout git branch `experiment/abc-smc-p0052-homeostasis-1000cells`.
- [ ] Add division inhibition contract methods to `ICell.h`.
- [ ] Implement `division_inhibited_` in `AgenticCell.h` and gate `phase5_Cytokinesis()`.
- [ ] Add `max_healthy_cells_` field and getter/setter in `Tissue.h`.
- [ ] Implement contact inhibition and daughter filtering in `Tissue::live()`.
- [ ] Add CLI flag `--max-healthy-cells` to `cell_scale/app/run_bootstrapping.cpp`.
- [ ] Create `cell_scale/tests/TissueHomeostasisTest.cpp` testing carrying capacity and tumoral escape.
- [ ] Verify the changes in terms of typechecking, linting and tests using the project's verification command. Fix issues if any.
- [ ] STOP. Present the changes to the user for review and suggest commit messages. Do NOT proceed to the next phase until the user explicitly asks.

---

### Phase 2: ABC-SMC Pipeline with Coupled Division and Broadened Priors

#### Description
Create a dedicated experiment directory `experiments/exp08_cs_abc_smc_homeostasis/` containing the revised ABC-SMC calibration runner. Incorporate basal division rate sampling $p \in [0.047, 0.057]$, strictly bind neoplastic division rate to $\nu = 10 \times p$, widen somatic mutation priors, configure 1000 starting cells, and evaluate against Kuchenbaecker et al. (2017) clinical cumulative incidence points.

#### Public Contracts
- **`experiments/exp08_cs_abc_smc_homeostasis/abc_smc.py`**:
  - `PARAM_NAMES = ["division_rate", "brca1_rate", "low_delta", "high_delta"]` (plus optional `tp53_rate` if varied).
  - Priors:
    - `division_rate`: Uniform $\mathcal{U}[0.047, 0.057]$
    - `brca1_rate`: Broadened uniform prior $\mathcal{U}[0.001, 0.150]$
    - `low_delta`: Uniform $\mathcal{U}[0.060, 0.220]$
    - `high_delta`: Uniform $\mathcal{U}[0.120, 0.600]$
  - Coupling: `--neoplastic-div-rate` passed dynamically as `10.0 * p_dict["division_rate"]`.
  - Fixed flags passed to `run_bootstrapping`:
    - `--cells 1000`
    - `--max-healthy-cells 1000`
    - `--big-bang`
    - `--weighted-sse`
  - Output files:
    - `generation_{gen:02d}.csv`
    - `smc_diagnostics.csv`

#### To-do Actions List
- [ ] Scaffold `experiments/exp08_cs_abc_smc_homeostasis/` directory.
- [ ] Implement `abc_smc.py` with updated parameter vector, priors, and dynamic $\nu = 10 \times p$ flag injection.
- [ ] Ensure generation CSV serialization records both `division_rate` and the coupled `neoplastic_div_rate`.
- [ ] Add CLI flags to `abc_smc.py` for configuring population size (default 1000), seeds, and generations.
- [ ] Verify the changes in terms of typechecking, linting and tests using the project's verification command. Fix issues if any.
- [ ] STOP. Present the changes to the user for review and suggest commit messages. Do NOT proceed to the next phase until the user explicitly asks.

---

### Phase 3: Pilot Validation Run and Diagnostics

#### Description
Execute a reduced proof-of-concept run of the new ABC-SMC pipeline (e.g. 2–3 generations, small particle count) to verify numerical stability, proper acceptance rates, and homeostatic maintenance of the healthy cell population. Document the findings and prepare the full run command.

#### Public Contracts
- **`experiments/exp08_cs_abc_smc_homeostasis/run_pilot.sh`**:
  - Script to run a fast proof-of-concept calibration.
- **`experiments/exp08_cs_abc_smc_homeostasis/README.md`**:
  - Documentation of the experiment design, rationale from `brca1_onset_dinamica_tejido.pdf`, and execution guide.

#### To-do Actions List
- [ ] Create `run_pilot.sh` to trigger a fast test run (e.g. 20 particles, 3 generations).
- [ ] Execute the pilot run and inspect `smc_diagnostics.csv` to ensure acceptable distance reduction without parameter collapse.
- [ ] Verify healthy cell counts stay bounded around 1000 throughout simulated tissue lifespans.
- [ ] Author `README.md` detailing the biological rationale, coupling rule, and instructions for full calibration.
- [ ] Verify the changes in terms of typechecking, linting and tests using the project's verification command. Fix issues if any.
- [ ] STOP. Present the changes to the user for review and suggest commit messages. Do NOT proceed to the next phase until the user explicitly asks.

---

## Next Step
Proceed to **Phase 1: Git Branch Creation and C++ Domain Homeostasis** upon user confirmation.
