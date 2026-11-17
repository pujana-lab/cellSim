#!/usr/bin/env python3
"""
plot_cell_stage_dynamics.py — Publication-grade visualization of cellular population dynamics
and mutation stage distributions over time (ages 0 to 80).

Panels:
  Panel a: Cellular Life Stages over time (Stacked area / Median curves: BASELINE, UNSTABLE, UNPROTECTED, PRIMER, TUMORAL).
  Panel b: Genotypic mutation dynamics (BRCA1 +/- vs -/-, TP53 +/+ vs +/- vs -/-).
  Panel c: Cell stage distribution at key clinical milestone ages (Ages 20, 40, 60, 80).

Outputs:
  - experiments/exp08_cs_abc_smc_homeostasis/Figure_stage_dynamics.pdf / .png
  - results/paper_figures/Figure_stage_dynamics.pdf / .png
"""

from __future__ import annotations

from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

plt.switch_backend("Agg")

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans", "Helvetica", "Arial"],
    "mathtext.fontset": "dejavusans",
    "axes.edgecolor": "#111827",
    "axes.linewidth": 0.8,
    "xtick.major.width": 0.8,
    "ytick.major.width": 0.8,
    "xtick.color": "#111827",
    "ytick.color": "#111827",
})


def plot_dynamics(csv_path: Path, out_dir: Path):
    df = pd.read_csv(csv_path)

    # Compute year-by-year mean and percentiles across tissues
    grouped = df.groupby("year")
    years = np.array(list(grouped.groups.keys()))

    stages = ["baseline", "unstable", "unprotected", "primer", "tumoral"]
    stage_labels = ["Baseline (TP53+/+, BRCA1+/-)", "Unstable (TP53+/-, BRCA1+/-)", 
                    "Unprotected (TP53-/-)", "Primer (D1 > 2.0)", "Tumoral (Transformed)"]
    stage_colors = ["#2563EB", "#F59E0B", "#10B981", "#8B5CF6", "#DC2626"]

    mean_stages = {s: grouped[s].mean().to_numpy() for s in stages}
    mean_total = grouped["total_alive"].mean().to_numpy()

    fig = plt.figure(figsize=(13.5, 9.0), dpi=300)
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 0.9], hspace=0.32, wspace=0.22)

    # =========================================================================
    # Panel a: Absolute cell counts by stage over time (Stacked Area)
    # =========================================================================
    ax_a = fig.add_subplot(gs[0, 0])
    fig.text(0.04, 0.96, "a", fontsize=15, fontweight="bold", va="top")
    ax_a.set_title("Evolution of Cellular Life Stages over Time (Mean Cell Count)", fontsize=11, fontweight="bold", pad=8)

    y_stack = [mean_stages[s] for s in stages]
    ax_a.stackplot(years, y_stack, labels=stage_labels, colors=stage_colors, alpha=0.88, edgecolor="#111827", linewidth=0.5)

    ax_a.plot(years, mean_total, color="#111827", linestyle="--", linewidth=1.5, label="Total Live Cells")
    ax_a.set_xlim(1, 80)
    ax_a.set_ylim(0, max(mean_total) * 1.15)
    ax_a.set_xlabel("Tissue Age (Years)", fontsize=10)
    ax_a.set_ylabel("Live Cells per Tissue", fontsize=10)
    ax_a.tick_params(labelsize=9)
    ax_a.legend(loc="upper right", fontsize=8.0, frameon=True, edgecolor="#E5E7EB")

    # =========================================================================
    # Panel b: Proportion (%) of Cell Stages over Time (Relative Composition)
    # =========================================================================
    ax_b = fig.add_subplot(gs[0, 1])
    fig.text(0.52, 0.96, "b", fontsize=15, fontweight="bold", va="top")
    ax_b.set_title("Relative Composition of Cellular Stages (%)", fontsize=11, fontweight="bold", pad=8)

    # Normalized proportions
    prop_stack = [mean_stages[s] / np.maximum(mean_total, 1e-9) * 100.0 for s in stages]
    ax_b.stackplot(years, prop_stack, labels=stage_labels, colors=stage_colors, alpha=0.88, edgecolor="#111827", linewidth=0.5)

    ax_b.set_xlim(1, 80)
    ax_b.set_ylim(0, 100)
    ax_b.set_xlabel("Tissue Age (Years)", fontsize=10)
    ax_b.set_ylabel("Percentage of Live Population (%)", fontsize=10)
    ax_b.tick_params(labelsize=9)
    ax_b.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), fontsize=8.0, frameon=False)

    # =========================================================================
    # Panel c: Genotype status dynamics (TP53 and BRCA1)
    # =========================================================================
    ax_c = fig.add_subplot(gs[1, 0])
    fig.text(0.04, 0.48, "c", fontsize=15, fontweight="bold", va="top")
    ax_c.set_title("Genetic Allele Dynamics (BRCA1 & TP53 Breakdown)", fontsize=11, fontweight="bold", pad=8)

    mean_b_pm = grouped["brca1_pm"].mean().to_numpy() / np.maximum(mean_total, 1e-9) * 100.0
    mean_b_mm = grouped["brca1_mm"].mean().to_numpy() / np.maximum(mean_total, 1e-9) * 100.0

    mean_tp_pp = grouped["tp53_pp"].mean().to_numpy() / np.maximum(mean_total, 1e-9) * 100.0
    mean_tp_pm = grouped["tp53_pm"].mean().to_numpy() / np.maximum(mean_total, 1e-9) * 100.0
    mean_tp_mm = grouped["tp53_mm"].mean().to_numpy() / np.maximum(mean_total, 1e-9) * 100.0

    ax_c.plot(years, mean_b_pm, color="#3B82F6", lw=1.8, label="BRCA1 +/- (Carrier)")
    ax_c.plot(years, mean_b_mm, color="#1D4ED8", lw=1.8, linestyle="--", label="BRCA1 -/- (LOH)")
    ax_c.plot(years, mean_tp_pp, color="#10B981", lw=1.8, label="TP53 +/+ (Wildtype)")
    ax_c.plot(years, mean_tp_pm, color="#F59E0B", lw=1.8, label="TP53 +/- (Heterozygous)")
    ax_c.plot(years, mean_tp_mm, color="#DC2626", lw=1.8, linestyle="--", label="TP53 -/- (Double Null)")

    ax_c.set_xlim(1, 80)
    ax_c.set_ylim(0, 105)
    ax_c.set_xlabel("Tissue Age (Years)", fontsize=10)
    ax_c.set_ylabel("Fraction of Total Cells (%)", fontsize=10)
    ax_c.tick_params(labelsize=9)
    ax_c.legend(loc="upper right", fontsize=8.0, frameon=True, edgecolor="#E5E7EB")

    # =========================================================================
    # Panel d: Distribution at Clinical Milestones (Ages 20, 40, 60, 80)
    # =========================================================================
    ax_d = fig.add_subplot(gs[1, 1])
    fig.text(0.52, 0.48, "d", fontsize=15, fontweight="bold", va="top")
    ax_d.set_title("Stage Distribution at Milestone Ages", fontsize=11, fontweight="bold", pad=8)

    milestone_ages = [10, 30, 50, 70]
    bar_width = 0.55
    x_positions = np.arange(len(milestone_ages))

    bottom_vals = np.zeros(len(milestone_ages))
    for s_idx, s in enumerate(stages):
        vals = []
        for age in milestone_ages:
            sub = df[df["year"] == age]
            tot = sub["total_alive"].mean()
            val = (sub[s].mean() / max(tot, 1e-9)) * 100.0
            vals.append(val)
        vals = np.array(vals)
        ax_d.bar(x_positions, vals, bar_width, bottom=bottom_vals, label=stage_labels[s_idx],
                 color=stage_colors[s_idx], edgecolor="#111827", lw=0.5)
        bottom_vals += vals

    ax_d.set_xticks(x_positions)
    ax_d.set_xticklabels([f"Age {age}" for age in milestone_ages], fontsize=9.5)
    ax_d.set_ylim(0, 100)
    ax_d.set_xlabel("Clinical Milestone Ages", fontsize=10)
    ax_d.set_ylabel("Population Fraction (%)", fontsize=10)
    ax_d.tick_params(labelsize=9)

    out_dir.mkdir(parents=True, exist_ok=True)
    pdf_path = out_dir / "Figure_stage_dynamics.pdf"
    png_path = out_dir / "Figure_stage_dynamics.png"

    fig.savefig(pdf_path, bbox_inches="tight")
    fig.savefig(png_path, dpi=300, bbox_inches="tight")
    plt.close(fig)

    # Mirror to results/paper_figures
    paper_dir = Path("results/paper_figures")
    paper_dir.mkdir(parents=True, exist_ok=True)
    import shutil
    shutil.copy2(pdf_path, paper_dir / "Figure_stage_dynamics.pdf")
    shutil.copy2(png_path, paper_dir / "Figure_stage_dynamics.png")

    print(f"Saved: {pdf_path}")
    print(f"Saved: {png_path}")
    print(f"Mirrored to: {paper_dir / 'Figure_stage_dynamics.png'}")


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Plot cell stage dynamics")
    parser.add_argument("--csv", type=str, default="experiments/exp08_cs_abc_smc_homeostasis/data/stage_trajectories_100tissues.csv", help="Input trajectories CSV")
    parser.add_argument("--out-dir", type=str, default="experiments/exp08_cs_abc_smc_homeostasis", help="Output directory")
    args = parser.parse_args()

    csv_file = Path(args.csv)
    out = Path(args.out_dir)
    plot_dynamics(csv_file, out)
