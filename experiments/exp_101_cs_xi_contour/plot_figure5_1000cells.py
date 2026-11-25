#!/usr/bin/env python3
"""
plot_figure5_1000cells.py — Publication-grade visualization for Figure 5 and Supplementary Figure S5.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Tuple

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

COLORS = {
    "none": "#1F77B4",
    "D1": "#FF7F0E",
    "D2": "#2CA02C",
    "D1+D2": "#D62728",
}

LABELS = {
    "none": "No intervention",
    "D1": r"$\theta_{D1}$ attenuation",
    "D2": r"$\theta_{D2}$ attenuation",
    "D1+D2": r"$\theta_{D1} + \theta_{D2}$ attenuation",
}


def render_figure5(data_dir: Path, out_dir: Path) -> Tuple[Path, Path]:
    df_curves = pd.read_csv(data_dir / "Figure5a_risk_curves.csv")
    df_arr = pd.read_csv(data_dir / "Figure5b_ARR_age70_50pct.csv")
    df_dose = pd.read_csv(data_dir / "Figure5c_dose_response_age30.csv")

    fig = plt.figure(figsize=(11.0, 4.2), dpi=300)
    gs = fig.add_gridspec(1, 3, left=0.06, right=0.96, top=0.84, bottom=0.18, wspace=0.35)

    fig.text(0.04, 0.96, "Figure 5", fontsize=14, fontweight="bold", va="top")
    fig.text(0.12, 0.96, "Cancer Interception (1000 cells neutral dynamics, F = 1.0)", fontsize=11, color="#374151", va="top")

    # Panel a
    ax_a = fig.add_subplot(gs[0, 0])
    ax_a.text(-0.20, 1.12, "a", transform=ax_a.transAxes, fontsize=13, fontweight="bold", va="top")
    ax_a.set_title("Modeled interception from age 30", fontsize=10.5, pad=6)

    for cond in ["none", "D1", "D2", "D1+D2"]:
        sub = df_curves[df_curves["condition"] == cond]
        ax_a.plot(sub["age"], sub["risk"], color=COLORS[cond], label=LABELS[cond], linewidth=1.4)

    ax_a.set_xlim(30, 80)
    ax_a.set_ylim(0, 80)
    ax_a.set_xticks([30, 40, 50, 60, 70, 80])
    ax_a.set_yticks([0, 20, 40, 60, 80])
    ax_a.set_xlabel("Age (years)", fontsize=9.5)
    ax_a.set_ylabel("Subsequent tumor-positive tissues (%)", fontsize=9.5)
    ax_a.tick_params(labelsize=8.5, length=2.5)
    ax_a.grid(True, linestyle=":", alpha=0.5, color="#E5E7EB")
    ax_a.legend(loc="upper left", fontsize=7.5, frameon=True, edgecolor="#E5E7EB", framealpha=0.9)

    # Panel b
    ax_b = fig.add_subplot(gs[0, 1])
    ax_b.text(-0.20, 1.12, "b", transform=ax_b.transAxes, fontsize=13, fontweight="bold", va="top")
    ax_b.set_title("Effect of intervention timing", fontsize=10.5, pad=6)

    start_ages = [30, 40, 50]
    x = np.arange(len(start_ages))
    width = 0.24

    for idx, cond in enumerate(["D1", "D2", "D1+D2"]):
        sub = df_arr[df_arr["condition"] == cond]
        arr_vals = [sub[sub["intervention_age"] == a]["ARR_pp"].values[0] for a in start_ages]
        ax_b.bar(x + (idx - 1) * width, arr_vals, width, label=LABELS[cond],
                 color=COLORS[cond], edgecolor="#111827", linewidth=0.5)

    ax_b.set_xticks(x)
    ax_b.set_xticklabels([f"Age {a}" for a in start_ages], fontsize=8.5)
    ax_b.set_xlabel("Intervention start age", fontsize=9.5)
    ax_b.set_ylabel("ARR by age 70 (% points)", fontsize=9.5)
    max_arr = df_arr["ARR_pp"].max()
    ax_b.set_ylim(0, max_arr * 1.25 if max_arr > 0 else 50)
    ax_b.tick_params(labelsize=8.5, length=2.5)
    ax_b.grid(True, linestyle=":", alpha=0.5, color="#E5E7EB", axis="y")
    ax_b.legend(loc="upper right", fontsize=7.0, frameon=True, edgecolor="#E5E7EB", framealpha=0.9)

    # Panel c
    ax_c = fig.add_subplot(gs[0, 2])
    ax_c.text(-0.20, 1.12, "c", transform=ax_c.transAxes, fontsize=13, fontweight="bold", va="top")
    ax_c.set_title("Dose-dependent interception effects", fontsize=10.5, pad=6)

    doses = [0.10, 0.25, 0.50]
    xd = np.arange(len(doses))

    for idx, cond in enumerate(["D1", "D2", "D1+D2"]):
        sub = df_dose[df_dose["condition"] == cond]
        arr_vals = [sub[np.isclose(sub["dose"], d)]["ARR_pp"].values[0] for d in doses]
        ax_c.bar(xd + (idx - 1) * width, arr_vals, width, label=LABELS[cond],
                 color=COLORS[cond], edgecolor="#111827", linewidth=0.5)

    ax_c.set_xticks(xd)
    ax_c.set_xticklabels(["10%", "25%", "50%"], fontsize=8.5)
    ax_c.set_xlabel("Intervention dose (age 30 start)", fontsize=9.5)
    ax_c.set_ylabel("ARR by age 70 (% points)", fontsize=9.5)
    max_dose_arr = df_dose["ARR_pp"].max()
    ax_c.set_ylim(0, max_dose_arr * 1.25 if max_dose_arr > 0 else 50)
    ax_c.tick_params(labelsize=8.5, length=2.5)
    ax_c.grid(True, linestyle=":", alpha=0.5, color="#E5E7EB", axis="y")
    ax_c.legend(loc="upper left", fontsize=7.0, frameon=True, edgecolor="#E5E7EB", framealpha=0.9)

    out_png = out_dir / "Figure5.png"
    out_pdf = out_dir / "Figure5.pdf"
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.savefig(out_pdf, bbox_inches="tight")
    plt.close()

    print(f"Figure 5 guardada en:\n  -> {out_png}")
    return out_pdf, out_png


def render_figures5(data_dir: Path, out_dir: Path):
    df_inter = pd.read_csv(data_dir / "FigureS5_interaction.csv")

    fig, ax = plt.subplots(figsize=(6.5, 4.2), dpi=300)
    fig.text(0.04, 0.96, "Supplementary Figure S5", fontsize=13, fontweight="bold", va="top")
    ax.set_title(r"Interaction Contrast on Absolute Risk ($D_1 + D_2$ vs Independent Additivity)", fontsize=10, pad=8)

    start_ages = [30, 40, 50]
    doses = [0.25, 0.50]
    x = np.arange(len(start_ages))
    width = 0.35

    dose_colors = {0.25: "#3B82F6", 0.50: "#1D4ED8"}

    for idx, dose in enumerate(doses):
        sub = df_inter[np.isclose(df_inter["dose"], dose)]
        vals = [sub[sub["intervention_age"] == a]["interaction_pp"].values[0] for a in start_ages]
        ci_lo = [sub[sub["intervention_age"] == a]["ci_lo"].values[0] for a in start_ages]
        ci_hi = [sub[sub["intervention_age"] == a]["ci_hi"].values[0] for a in start_ages]
        yerr = [
            [v - l for v, l in zip(vals, ci_lo)],
            [h - v for v, h in zip(vals, ci_hi)]
        ]
        ax.bar(x + (idx - 0.5) * width, vals, width, yerr=yerr, capsize=4,
               label=f"{int(dose*100)}% attenuation", color=dose_colors[dose], edgecolor="#111827", linewidth=0.6)

    ax.axhline(0, color="#DC2626", linestyle="--", linewidth=1.0, alpha=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels([f"Age {a}" for a in start_ages], fontsize=9)
    ax.set_xlabel("Intervention start age", fontsize=9.5)
    ax.set_ylabel("Interaction contrast (% points through age 80)", fontsize=9.5)
    ax.tick_params(labelsize=8.5)
    ax.grid(True, linestyle=":", alpha=0.5, color="#E5E7EB", axis="y")
    ax.legend(loc="lower left", fontsize=8, frameon=True, edgecolor="#E5E7EB")

    out_png = out_dir / "FigureS5.png"
    out_pdf = out_dir / "FigureS5.pdf"
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.savefig(out_pdf, bbox_inches="tight")
    plt.close()
    print(f"Figure S5 guardada en:\n  -> {out_png}")


def render_figure5_4panels(data_dir: Path, out_dir: Path) -> Tuple[Path, Path]:
    df_curves_50 = pd.read_csv(data_dir / "Figure5a_risk_curves_50pct.csv")
    df_curves_75 = pd.read_csv(data_dir / "Figure5b_risk_curves_75pct.csv")
    df_arr_75 = pd.read_csv(data_dir / "Figure5c_ARR_age70_75pct.csv")
    df_dose = pd.read_csv(data_dir / "Figure5d_dose_response_age30.csv")

    fig = plt.figure(figsize=(15.0, 4.2), dpi=300)
    gs = fig.add_gridspec(1, 4, left=0.05, right=0.98, top=0.84, bottom=0.18, wspace=0.35)

    fig.text(0.03, 0.96, "Figure 5", fontsize=14, fontweight="bold", va="top")
    fig.text(0.10, 0.96, "Cancer Interception (1000 cells neutral dynamics, F = 1.0)", fontsize=11, color="#374151", va="top")

    # Panel a: 50% attenuation (2.0x fold-reduction)
    ax_a = fig.add_subplot(gs[0, 0])
    ax_a.text(-0.20, 1.12, "a", transform=ax_a.transAxes, fontsize=13, fontweight="bold", va="top")
    ax_a.set_title("Interception age 30 (50% dose, 2.0× FC)", fontsize=9.5, pad=6)
    for cond in ["none", "D1", "D2", "D1+D2"]:
        sub = df_curves_50[df_curves_50["condition"] == cond]
        ax_a.plot(sub["age"], sub["risk"], color=COLORS[cond], label=LABELS[cond], linewidth=1.4)
    ax_a.set_xlim(30, 80)
    ax_a.set_ylim(0, 80)
    ax_a.set_xticks([30, 40, 50, 60, 70, 80])
    ax_a.set_yticks([0, 20, 40, 60, 80])
    ax_a.set_xlabel("Age (years)", fontsize=9.0)
    ax_a.set_ylabel("Subsequent tumor-positive tissues (%)", fontsize=9.0)
    ax_a.tick_params(labelsize=8.0, length=2.5)
    ax_a.grid(True, linestyle=":", alpha=0.5, color="#E5E7EB")
    ax_a.legend(loc="upper left", fontsize=6.8, frameon=True, edgecolor="#E5E7EB", framealpha=0.9)

    # Panel b: 75% attenuation (4.0x fold-reduction)
    ax_b = fig.add_subplot(gs[0, 1])
    ax_b.text(-0.20, 1.12, "b", transform=ax_b.transAxes, fontsize=13, fontweight="bold", va="top")
    ax_b.set_title("Interception age 30 (75% dose, 4.0× FC)", fontsize=9.5, pad=6)
    for cond in ["none", "D1", "D2", "D1+D2"]:
        sub = df_curves_75[df_curves_75["condition"] == cond]
        ax_b.plot(sub["age"], sub["risk"], color=COLORS[cond], label=LABELS[cond], linewidth=1.4)
    ax_b.set_xlim(30, 80)
    ax_b.set_ylim(0, 80)
    ax_b.set_xticks([30, 40, 50, 60, 70, 80])
    ax_b.set_yticks([0, 20, 40, 60, 80])
    ax_b.set_xlabel("Age (years)", fontsize=9.0)
    ax_b.set_ylabel("Subsequent tumor-positive tissues (%)", fontsize=9.0)
    ax_b.tick_params(labelsize=8.0, length=2.5)
    ax_b.grid(True, linestyle=":", alpha=0.5, color="#E5E7EB")
    ax_b.legend(loc="upper left", fontsize=6.8, frameon=True, edgecolor="#E5E7EB", framealpha=0.9)

    # Panel c: Effect of intervention timing at 75% dose (4.0x fold-reduction)
    ax_c = fig.add_subplot(gs[0, 2])
    ax_c.text(-0.20, 1.12, "c", transform=ax_c.transAxes, fontsize=13, fontweight="bold", va="top")
    ax_c.set_title("Timing effect (75% dose, 4.0× FC)", fontsize=9.5, pad=6)
    start_ages = [30, 40, 50]
    x = np.arange(len(start_ages))
    width = 0.24
    for idx, cond in enumerate(["D1", "D2", "D1+D2"]):
        sub = df_arr_75[df_arr_75["condition"] == cond]
        arr_vals = [sub[sub["intervention_age"] == a]["ARR_pp"].values[0] for a in start_ages]
        ax_c.bar(x + (idx - 1) * width, arr_vals, width, label=LABELS[cond],
                 color=COLORS[cond], edgecolor="#111827", linewidth=0.5)
    ax_c.set_xticks(x)
    ax_c.set_xticklabels([f"Age {a}" for a in start_ages], fontsize=8.0)
    ax_c.set_xlabel("Intervention start age", fontsize=9.0)
    ax_c.set_ylabel("ARR by age 70 (% points)", fontsize=9.0)
    max_arr = df_arr_75["ARR_pp"].max()
    ax_c.set_ylim(0, max_arr * 1.25 if max_arr > 0 else 50)
    ax_c.tick_params(labelsize=8.0, length=2.5)
    ax_c.grid(True, linestyle=":", alpha=0.5, color="#E5E7EB", axis="y")
    ax_c.legend(loc="upper right", fontsize=6.8, frameon=True, edgecolor="#E5E7EB", framealpha=0.9)

    # Panel d: Dose-dependent interception effects (10, 25, 50, 75%) with Fold-Change
    ax_d = fig.add_subplot(gs[0, 3])
    ax_d.text(-0.20, 1.12, "d", transform=ax_d.transAxes, fontsize=13, fontweight="bold", va="top")
    ax_d.set_title("Dose-dependent effects (age 30 start)", fontsize=9.5, pad=6)
    doses = [0.10, 0.25, 0.50, 0.75]
    xd = np.arange(len(doses))
    for idx, cond in enumerate(["D1", "D2", "D1+D2"]):
        sub = df_dose[df_dose["condition"] == cond]
        arr_vals = [sub[np.isclose(sub["dose"], d)]["ARR_pp"].values[0] for d in doses]
        ax_d.bar(xd + (idx - 1) * width, arr_vals, width, label=LABELS[cond],
                 color=COLORS[cond], edgecolor="#111827", linewidth=0.5)
    ax_d.set_xticks(xd)
    ax_d.set_xticklabels([
        "10%\n(1.1×)",
        "25%\n(1.3×)",
        "50%\n(2.0×)",
        "75%\n(4.0×)"
    ], fontsize=8.0)
    ax_d.set_xlabel("Intervention dose [Fold-reduction]", fontsize=9.0)
    ax_d.set_ylabel("ARR by age 70 (% points)", fontsize=9.0)
    max_dose_arr = df_dose["ARR_pp"].max()
    ax_d.set_ylim(0, max_dose_arr * 1.25 if max_dose_arr > 0 else 50)
    ax_d.tick_params(labelsize=8.0, length=2.5)
    ax_d.grid(True, linestyle=":", alpha=0.5, color="#E5E7EB", axis="y")
    ax_d.legend(loc="upper left", fontsize=6.8, frameon=True, edgecolor="#E5E7EB", framealpha=0.9)

    out_png = out_dir / "Figure5.png"
    out_pdf = out_dir / "Figure5.pdf"
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.savefig(out_pdf, bbox_inches="tight")
    plt.close()

    print(f"Figure 5 (4 paneles) guardada en:\n  -> {out_png}")
    return out_pdf, out_png


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", type=Path, default=Path("experiments/exp_101_cs_xi_contour/results_1000cells_final/interception"))
    parser.add_argument("--out-dir", type=Path, default=Path("experiments/exp_101_cs_xi_contour/results_1000cells_final"))
    args = parser.parse_args()

    render_figure5_4panels(args.data_dir, args.out_dir)
    render_figures5(args.data_dir, args.out_dir)


if __name__ == "__main__":
    main()
