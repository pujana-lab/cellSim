#!/usr/bin/env python3
"""
plot_figure4_homeostasis.py — Publication-grade visualization for Figure 4 using the
calibrated homeostatic ABC-SMC posterior (1000 cells, p in [0.047, 0.057], nu = 10*p).

Panels:
  Panel a: Calibrated parameter posteriors (division_rate, brca1_rate, low_delta, high_delta)
           - Row 1: 1D KDEs with uniform priors and medians.
           - Row 2: Empirical posterior CDFs vs theoretical uniform CDFs.
           - Right: Weighted pairwise correlation matrix.
  Panel b: Somatic BRCA1 status in emerging tumors stratified by detection age
           (BRCA1-/- null, Mixed, BRCA1+/- heterozygous).
  Panel c: Tumoral founder-lineage complexity by detection age (1, 2, 3, >=4 founders).
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Tuple

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

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

PARAMS = ["division_rate", "brca1_rate", "low_delta", "high_delta"]
PARAM_LABELS = {
    "division_rate": r"$p$",
    "brca1_rate": r"$\beta$",
    "low_delta": r"$\delta_{\rm low}$",
    "high_delta": r"$\delta_{\rm high}$",
}
PRIORS = {
    "division_rate": (0.047, 0.057),
    "brca1_rate": (0.001, 0.150),
    "low_delta": (0.060, 0.220),
    "high_delta": (0.120, 0.600),
}
TICKS = {
    "division_rate": [0.048, 0.052, 0.056],
    "brca1_rate": [0.02, 0.08, 0.14],
    "low_delta": [0.08, 0.14, 0.20],
    "high_delta": [0.2, 0.4, 0.6],
}


def weighted_stats(x: np.ndarray, w: np.ndarray) -> Tuple[float, float, float]:
    idx = np.argsort(x)
    sx = x[idx]
    sw = w[idx]
    cw = np.cumsum(sw) / np.sum(sw)
    q025 = float(np.interp(0.025, cw, sx))
    q50 = float(np.interp(0.50, cw, sx))
    q975 = float(np.interp(0.975, cw, sx))
    return q50, q025, q975


def render_figure4(post_file: Path, features_dir: Path, out_dir: Path) -> Tuple[Path, Path]:
    df_part = pd.read_csv(post_file)
    w = df_part["weight"].to_numpy(float, copy=True)
    w = w / np.sum(w)

    df_geno_age = pd.read_csv(features_dir / "Figure4b_genotype_composition_by_age.csv")
    df_found = pd.read_csv(features_dir / "Figure4c_founder_lineages_by_age.csv")

    fig = plt.figure(figsize=(10.8, 9.8), dpi=300)
    fig.text(0.04, 0.98, "Figure 4", fontsize=15, fontweight="bold", va="top")

    # =========================================================================
    # Panel a: Calibrated parameter posteriors
    # =========================================================================
    fig.text(0.04, 0.94, "a", fontsize=14, fontweight="bold", va="top")
    fig.text(0.07, 0.94, "Calibrated parameter posteriors (1000 cells neutral dynamics)", fontsize=11, fontweight="bold", va="top")

    gs_a = fig.add_gridspec(
        nrows=2, ncols=5,
        left=0.07, right=0.96, top=0.90, bottom=0.56,
        width_ratios=[1.0, 1.0, 1.0, 1.0, 1.15],
        wspace=0.38, hspace=0.35
    )

    y_max_dens = [250, 90, 40, 4.0]

    for col, p in enumerate(PARAMS):
        vals = df_part[p].to_numpy(float)
        p_lo, p_hi = PRIORS[p]
        grid_x = np.linspace(p_lo, p_hi, 200)

        # Row 0: Density KDE
        ax_dens = fig.add_subplot(gs_a[0, col])
        kde = stats.gaussian_kde(vals, weights=w)
        dens = kde(grid_x)
        ax_dens.plot(grid_x, dens, color="#111827", linewidth=1.3, label="Posterior KDE")

        prior_h = 1.0 / (p_hi - p_lo)
        ax_dens.axhline(prior_h, color="#9CA3AF", linestyle="--", linewidth=1.1, label="Uniform prior")

        med, _, _ = weighted_stats(vals, w)
        ax_dens.axvline(med, color="#DC2626", linewidth=1.3, label="Median")
        ax_dens.text(med, y_max_dens[col] * 0.90, f"{med:.3f}", color="#DC2626",
                     fontsize=8, ha="center", fontweight="bold")

        ax_dens.set_title(PARAM_LABELS[p], fontsize=12, pad=4)
        ax_dens.set_xlim(p_lo, p_hi)
        ax_dens.set_ylim(0, y_max_dens[col])
        ax_dens.set_xticks(TICKS[p])
        ax_dens.tick_params(labelsize=8, length=2.5)
        ax_dens.set_ylabel("Density" if col == 0 else "", fontsize=9)
        if col == 3:
            ax_dens.legend(loc="upper right", bbox_to_anchor=(1.35, 1.0), fontsize=7.5, frameon=False)

        # Row 1: CDF
        ax_cdf = fig.add_subplot(gs_a[1, col])
        idx_s = np.argsort(vals)
        s_vals = vals[idx_s]
        s_w = w[idx_s]
        ecdf = np.cumsum(s_w)
        ax_cdf.plot(s_vals, ecdf, color="#111827", linewidth=1.3, label="Posterior CDF")

        u_cdf = (grid_x - p_lo) / (p_hi - p_lo)
        ax_cdf.plot(grid_x, u_cdf, color="#9CA3AF", linestyle="--", linewidth=1.1, label="Theoretical prior CDF")

        ax_cdf.set_xlim(p_lo, p_hi)
        ax_cdf.set_ylim(0.0, 1.05)
        ax_cdf.set_xticks(TICKS[p])
        ax_cdf.set_yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
        ax_cdf.tick_params(labelsize=8, length=2.5)
        ax_cdf.set_ylabel("CDF" if col == 0 else "", fontsize=9)
        if col == 3:
            ax_cdf.legend(loc="upper right", bbox_to_anchor=(1.60, 1.0), fontsize=7.5, frameon=False)

    # Weighted Correlation Matrix
    ax_corr = fig.add_subplot(gs_a[:, 4])
    X = df_part[PARAMS].to_numpy(float)
    mean_w = np.average(X, weights=w, axis=0)
    X_c = X - mean_w
    cov_w = np.dot(w * X_c.T, X_c)
    diag_std = np.sqrt(np.diag(cov_w))
    corr_mat = cov_w / np.outer(diag_std, diag_std)

    im = ax_corr.imshow(corr_mat, cmap="coolwarm", vmin=-1.0, vmax=1.0)
    cbar = fig.colorbar(im, ax=ax_corr, fraction=0.046, pad=0.06)
    cbar.ax.tick_params(labelsize=8)
    cbar.set_ticks([-1.0, 0.0, 1.0])
    cbar.set_ticklabels(["-1", "0", "+1"])
    ax_corr.set_title("Weighted\ncorrelation", fontsize=9, pad=8)
    labels = [r"$p$", r"$\beta$", r"$\delta_{\rm low}$", r"$\delta_{\rm high}$"]
    ax_corr.set_xticks(range(4))
    ax_corr.set_yticks(range(4))
    ax_corr.set_xticklabels(labels, fontsize=10)
    ax_corr.set_yticklabels(labels, fontsize=10)
    ax_corr.tick_params(length=2.0)

    # =========================================================================
    # Panel b: BRCA1 status in emerging tumors (age-stratified)
    # =========================================================================
    fig.text(0.04, 0.48, "b", fontsize=14, fontweight="bold", va="top")
    fig.text(0.07, 0.48, "BRCA1 loss status in emerging tumors", fontsize=11, fontweight="bold", va="top")

    ax_b = fig.add_axes([0.08, 0.12, 0.39, 0.30])
    colors_b = ["#1E3A8A", "#F59E0B", "#93C5FD"]

    ages_b = df_geno_age["age_group"].tolist()
    xb = np.arange(len(ages_b))
    width = 0.65
    pmm = df_geno_age["exclusively_brca1_mm_pct"].to_numpy(float)
    pmix = df_geno_age["mixed_pct"].to_numpy(float)
    ppm = df_geno_age["exclusively_brca1_pm_pct"].to_numpy(float)

    ax_b.bar(xb, pmm, width, color=colors_b[0], label="BRCA1-/- (null)", edgecolor="#111827", lw=0.5)
    ax_b.bar(xb, pmix, width, bottom=pmm, color=colors_b[1], label="Mixed", edgecolor="#111827", lw=0.5)
    ax_b.bar(xb, ppm, width, bottom=pmm + pmix, color=colors_b[2], label="BRCA1+/- (heterozygous)", edgecolor="#111827", lw=0.5)

    ax_b.set_xticks(xb)
    ax_b.set_xticklabels(["<30", "30–39", "40–49", "50–59", "60–69", "70–79", "80"], fontsize=8.5)
    ax_b.set_xlabel("Age at first tumor-positive detection (years)", fontsize=9.5)
    ax_b.set_ylabel("Tumor-positive tissues (%)", fontsize=9.5)
    ax_b.set_ylim(0, 105)
    ax_b.tick_params(length=2.5, labelsize=8.5)
    ax_b.legend(loc="upper right", fontsize=7.5, frameon=True, edgecolor="#E5E7EB")

    # =========================================================================
    # Panel c: Founder-lineage complexity by detection age
    # =========================================================================
    fig.text(0.52, 0.48, "c", fontsize=14, fontweight="bold", va="top")
    fig.text(0.55, 0.48, "Founder-lineage complexity by detection age", fontsize=11, fontweight="bold", va="top")

    ax_c = fig.add_axes([0.56, 0.12, 0.39, 0.30])
    ages = df_found["age_group"].tolist()
    p1 = df_found["founder_1_pct"].to_numpy(float)
    p2 = df_found["founder_2_pct"].to_numpy(float)
    p3 = df_found["founder_3_pct"].to_numpy(float)
    p4 = df_found["founder_>=4_pct"].to_numpy(float)

    x = np.arange(len(ages))
    colors_c = ["#F59E0B", "#F97316", "#DC2626", "#991B1B"]

    ax_c.bar(x, p1, width, color=colors_c[0], label="1")
    ax_c.bar(x, p2, width, bottom=p1, color=colors_c[1], label="2")
    ax_c.bar(x, p3, width, bottom=p1 + p2, color=colors_c[2], label="3")
    ax_c.bar(x, p4, width, bottom=p1 + p2 + p3, color=colors_c[3], label=r"$\geq 4$")

    ax_c.set_xticks(x)
    ax_c.set_xticklabels(["<30", "30–39", "40–49", "50–59", "60–69", "70–79", "80"], fontsize=8.5)
    ax_c.set_xlabel("Age at first tumor-positive detection (years)", fontsize=9.5)
    ax_c.set_ylabel("Tumor-positive tissues (%)", fontsize=9.5)
    ax_c.set_ylim(0, 100)
    ax_c.tick_params(length=2.5, labelsize=8.5)

    stats_file = features_dir / "stats_validation.txt"
    if stats_file.exists():
        with open(stats_file) as f:
            lines = dict(l.strip().split(": ") for l in f if ": " in l)
        tau_val = float(lines.get("Kendall_tau", -0.065))
        p_val = float(lines.get("Kendall_p", 1e-13))
        if p_val < 1e-12:
            stat_text = rf"Kendall's $\tau_b = {tau_val:.2f}$ ($P < 10^{{-12}}$)"
        else:
            stat_text = rf"Kendall's $\tau_b = {tau_val:.2f}$ ($P = {p_val:.1e}$)"
        ax_c.text(0.70, 0.68, stat_text, transform=ax_c.transAxes,
                  fontsize=8.5, fontstyle="italic", ha="center")

    ax_c.legend(title="Tumoral founder lineages (n)", title_fontsize=8, loc="upper right",
                fontsize=8, frameon=True, edgecolor="#E5E7EB")

    # Output paths
    pdf_path = out_dir / "Figure4.pdf"
    png_path = out_dir / "Figure4.png"
    fig.savefig(pdf_path, bbox_inches="tight")
    fig.savefig(png_path, dpi=300, bbox_inches="tight")
    plt.close(fig)

    return pdf_path, png_path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--post-file", type=str, default="experiments/exp08_cs_abc_smc_homeostasis/data/generation_14.csv")
    parser.add_argument("--features-dir", type=str, default="experiments/exp08_cs_abc_smc_homeostasis/emergent_features")
    parser.add_argument("--out-dir", type=str, default="experiments/exp08_cs_abc_smc_homeostasis")
    args = parser.parse_args()

    pdf_p, png_p = render_figure4(Path(args.post_file), Path(args.features_dir), Path(args.out_dir))
    print(f"[SUCCESS] Generated: {pdf_p} and {png_p}")


if __name__ == "__main__":
    main()
