#!/usr/bin/env python3
"""
explore_results.py — Summary and visual diagnostics for exp08 ABC-SMC calibration.
"""

from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

SCRIPT_DIR = Path(__file__).parent
RESULTS_DIR = SCRIPT_DIR / "results_pilot"

# Kuchenbaecker et al. (2017) clinical cumulative breast cancer incidence (JAMA 2017)
CLINICAL_AGES = [30, 40, 50, 60, 70, 80]
CLINICAL_RISK = [4.0, 26.0, 46.0, 58.0, 65.0, 70.0]
CLINICAL_LOWER = [2.0, 22.0, 41.0, 52.0, 56.0, 60.0]
CLINICAL_UPPER = [7.0, 30.0, 52.0, 65.0, 73.0, 80.0]


def main():
    print(f"\n{'='*75}")
    print(" EXP08 ABC-SMC PILOT: RESULT EXPLORATION & DIAGNOSTICS")
    print(f"{'='*75}\n")

    diag_file = RESULTS_DIR / "smc_diagnostics.csv"
    if diag_file.exists():
        df_diag = pd.read_csv(diag_file)
        print("1. SMC Convergence Diagnostics:")
        print(df_diag[["generation", "epsilon", "acceptance_rate", "ess", "weighted_sse_median", "weighted_sse_min"]].to_string(index=False))
        print()

    # Load Gen 0 and latest Gen
    gen0_file = RESULTS_DIR / "generation_00.csv"
    gen_files = sorted(RESULTS_DIR.glob("generation_[0-9][0-9].csv"))
    if not gen_files:
        print("No generation files found.")
        return

    latest_file = gen_files[-1]
    df_gen0 = pd.read_csv(gen0_file) if gen0_file.exists() else None
    df_latest = pd.read_csv(latest_file)

    print(f"2. Parameter Posteriors ({latest_file.name}):")
    params = ["division_rate", "brca1_rate", "low_delta", "high_delta", "neoplastic_div_rate"]
    w = df_latest["weight"].values
    w_norm = w / np.sum(w)

    summary_rows = []
    for p in params:
        vals = df_latest[p].values
        mean = np.sum(vals * w_norm)
        # Weighted percentiles
        sort_idx = np.argsort(vals)
        cum_w = np.cumsum(w_norm[sort_idx])
        p025 = vals[sort_idx][np.searchsorted(cum_w, 0.025)]
        median = vals[sort_idx][np.searchsorted(cum_w, 0.50)]
        p975 = vals[sort_idx][np.searchsorted(cum_w, 0.975)]
        summary_rows.append({
            "Parameter": p,
            "Weighted Mean": f"{mean:.5f}",
            "Median": f"{median:.5f}",
            "95% Credible Interval": f"[{p025:.5f}, {p975:.5f}]"
        })
    print(pd.DataFrame(summary_rows).to_string(index=False))
    print()

    # Best particle (minimum SSE)
    best_idx = df_latest["distance_weighted_sse"].idxmin()
    best_row = df_latest.loc[best_idx]
    print(f"3. Best Particle in {latest_file.name} (Weighted SSE = {best_row['distance_weighted_sse']:.2f}):")
    for p in params:
        print(f"   {p}: {best_row[p]:.5f}")
    print()

    # Plot Visual Diagnostics
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    # Panel A: Convergence epsilon & median SSE
    if diag_file.exists():
        ax = axes[0]
        gens = df_diag["generation"].values
        ax.plot(gens, df_diag["weighted_sse_median"], "o-", color="#1f77b4", label="Median weighted SSE", lw=2)
        ax.plot(gens, df_diag["weighted_sse_min"], "s--", color="#2ca02c", label="Min weighted SSE", lw=2)
        ax.set_xlabel("SMC Generation")
        ax.set_ylabel("Weighted SSE")
        ax.set_title("A: Objective Distance Reduction")
        ax.grid(True, alpha=0.3)
        ax.legend()

    # Panel B: Cumulative penetrance curves vs Kuchenbaecker
    ax = axes[1]
    # Clinical points with 95% CI
    yerr = [
        np.array(CLINICAL_RISK) - np.array(CLINICAL_LOWER),
        np.array(CLINICAL_UPPER) - np.array(CLINICAL_RISK)
    ]
    ax.errorbar(CLINICAL_AGES, CLINICAL_RISK, yerr=yerr, fmt="o", color="black",
                label="Kuchenbaecker (2017) 95% CI", capsize=4, zorder=5)

    # Plot sample particles from latest generation
    for _, row in df_latest.iterrows():
        risks = [row[f"risk_{age}"] for age in CLINICAL_AGES]
        ax.plot(CLINICAL_AGES, risks, color="#1f77b4", alpha=0.25)

    # Plot best particle
    best_risks = [best_row[f"risk_{age}"] for age in CLINICAL_AGES]
    ax.plot(CLINICAL_AGES, best_risks, color="#d62728", lw=2.5, label=f"Best particle (SSE={best_row['distance_weighted_sse']:.1f})")

    ax.set_xlabel("Age (years)")
    ax.set_ylabel("Cumulative Breast Cancer Risk (%)")
    ax.set_title("B: Cumulative Risk vs. Clinical Benchmark")
    ax.set_ylim(-2, 105)
    ax.grid(True, alpha=0.3)
    ax.legend()

    # Panel C: Division vs Mutation rate scatter (d vs beta)
    ax = axes[2]
    if df_gen0 is not None:
        ax.scatter(df_gen0["division_rate"], df_gen0["brca1_rate"], color="gray", alpha=0.4, label="Prior (Gen 0)")
    ax.scatter(df_latest["division_rate"], df_latest["brca1_rate"], color="#d62728", alpha=0.8,
               s=df_latest["weight"] * 1500, label=f"Posterior ({latest_file.name})")
    ax.set_xlabel("Basal division rate p")
    ax.set_ylabel("BRCA1 mutation rate beta")
    ax.set_title("C: Division (p) vs Mutation (beta) Posterior")
    ax.grid(True, alpha=0.3)
    ax.legend()

    plt.tight_layout()
    output_png = SCRIPT_DIR / "pilot_diagnostics.png"
    plt.savefig(output_png, dpi=200)
    print(f"Visual diagnostics saved to: {output_png}\n")


if __name__ == "__main__":
    main()
