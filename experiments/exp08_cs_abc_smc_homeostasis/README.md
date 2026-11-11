# Experiment 08: ABC-SMC Calibration with Healthy Tissue Homeostasis & Proliferation Coupling

## Rationale
Based on the analytical derivations in `brca1_onset_dinamica_tejido.pdf`, the model is updated with:
1. **Healthy tissue carrying capacity**: Fixed at $N_{\text{normal}} = 1000$ cells with contact inhibition. Neoplastic clones proliferate autonomously under Big Bang mode.
2. **Basal division rate**: $p \in [0.047, 0.057]$.
3. **Tumoral division coupling**: $\nu = 10 \times p \in [0.47, 0.57]$.
4. **Broadened mutation priors**: $\beta \in [0.001, 0.150]$.

## Quick Start: Running the Full Calibration
To run the full 20-generation ABC-SMC calibration:
```bash
python3 experiments/exp08_cs_abc_smc_homeostasis/abc_smc.py \
    --particles 100 \
    --generations 15 \
    --n-sim 200 \
    --n-cells 1000 \
    --max-healthy-cells 1000 \
    --workers 8 \
    --output-dir experiments/exp08_cs_abc_smc_homeostasis/results
```

## Pilot Run Results (4 generations, 20 particles)
- Gen 0: Median weighted SSE = 1616.05
- Gen 1: Median weighted SSE = 820.69
- Gen 2: Median weighted SSE = 791.21
- Gen 3: Median weighted SSE = 519.09 (min = 169.50)
- Execution time: ~65 seconds across 4 generations.
