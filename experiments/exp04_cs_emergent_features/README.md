# Experiment 04: Emergent Tumor Features & Clonal Lineages (Figure 4)

## 🎯 Overview
This vertical slice replicates **Figure 4** of the BRCA1-ONSET manuscript (*Palomero et al., 2026*).

It evaluates emergent tumor properties recovered by the calibrated discrete agent-based model across the full posterior population ($250 \text{ particles} \times 80 \text{ virtual tissues} = 20,000 \text{ simulated tissues}$):
1. **Panel a**: Posterior parameter distributions, theoretical uniform priors, and weighted pairwise correlation matrix.
2. **Panel b**: Age-stratified genotype composition of emerging tumors ($BRCA1^{+/-}$ vs Mixed vs $BRCA1^{-/-}$).
3. **Panel c**: Heterogeneity in founder-lineage complexity across age at tumor detection (Kendall's $\tau_b = -0.25$, $P < 10^{-15}$; $\chi^2\ P = 3.3 \times 10^{-231}$).

---

## 🧬 Biological & Computational Results (Simulated $N=20,000$ Tissues)

### 1. Genotype Composition of Emerging Tumors by Age (Figure 4b)
- **Global Breakdown**:
  - **Exclusively $BRCA1^{-/-}$**: 58.65% (manuscript: ~60%)
  - **Mixed ($BRCA1^{+/-}$ and $BRCA1^{-/-}$)**: 28.29% (manuscript: ~29%)
  - **Exclusively $BRCA1^{+/-}$**: 13.07% (manuscript: ~11%)
  - **Total tissues with a $BRCA1^{-/-}$ component**: **86.94%** (matches the 80–90% frequency reported in clinical literature).
- **Age-Stratified Evolution**:
  - $<30$ years: 0.57% exclusively $BRCA1^{+/-}$, 41.98% mixed, 57.45% exclusively $BRCA1^{-/-}$
  - $40–49$ years: 7.80% exclusively $BRCA1^{+/-}$, 29.35% mixed, 62.85% exclusively $BRCA1^{-/-}$
  - $\ge 70$ years: $>62\%$ exclusively $BRCA1^{+/-}$ (predominantly $BRCA1^{+/-}$ maintenance with late somatic bypass).

### 2. Founder Clonal Complexity by Age (Figure 4c)
- Total tumor-positive tissues analyzed: **$N = 13,309$**.
- Monoclonal vs. Multiclonal dynamics:
  - $<30$ years: 60.7% single founder, 39.3% multiclonal ($\ge 2$ founders).
  - $30–39$ years: 63.6% single founder, 36.4% multiclonal.
  - $40–49$ years: 77.7% single founder, 22.3% multiclonal.
  - $50–59$ years: 88.5% single founder, 11.5% multiclonal.
  - $\ge 70$ years: $>98.3\%$ single founder (predominantly monoclonal).
- **Statistical Trend Evaluation**:
  - **Kendall's rank correlation (ordinal trend)**: $\tau_b = -0.2461$, **$P < 10^{-15}$** ($P = 4.23 \times 10^{-264}$).
  - **Contingency test (reference)**: $\chi^2 = 1,141.79$, $df = 18$, $P = 3.29 \times 10^{-231}$.

---

## 🚀 Execution Commands

```bash
# 1. Render Figure 4 (Instantaneous via simulated data)
make exp04-fig4

# 2. Fast Pipeline Smoke Test (< 0.1s)
make exp04-smoke

# 3. Full Cohort Simulation (Optional recomputation of 20,000 tissues, ~45s)
make exp04-simulate
```
