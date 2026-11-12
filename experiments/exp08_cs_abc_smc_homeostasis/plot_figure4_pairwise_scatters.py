#!/usr/bin/env python3
"""
plot_figure4_pairwise_scatters.py — Publication-grade pairwise scatter matrix
for calibrated ABC-SMC parameters under homeostasis.

Displays all pairwise scatter plots and bivariate density contours for:
  - p (division_rate)
  - beta (brca1_rate)
  - delta_low (low_delta)
  - delta_high (high_delta)

Outputs:
  - Figure4_pairwise_scatters.pdf
  - Figure4_pairwise_scatters.png
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
    "brca1_rate": [0.03, 0.08, 0.13],
    "low_delta": [0.08, 0.14, 0.20],
    "high_delta": [0.2, 0.4, 0.6],
}


def weighted_stats(x: np.ndarray, w: np.ndarray) -> Tuple[float, float, float]:
    idx = np.argsort(x)
    sx, sw = x[idx], w[idx]
    cw = np.cumsum(sw)
    q025 = sx[np.searchsorted(cw, 0.025)]
    med = sx[np.searchsorted(cw, 0.50)]
    q975 = sx[np.searchsorted(cw, 0.975)]
    return med, q025, q975


def render_pairwise_scatters(data_file: Path, out_dir: Path) -> Tuple[Path, Path]:
    if not data_file.exists():
        raise FileNotFoundError(f"Posterior data file not found: {data_file}")

    df = pd.read_csv(data_file)
    w = df["weight"].to_numpy(float)
    w = w / np.sum(w)

    n_params = len(PARAMS)
    fig, axes = plt.subplots(n_params, n_params, figsize=(8.5, 8.5), dpi=300)
    plt.subplots_adjust(left=0.10, right=0.96, bottom=0.09, top=0.95, hspace=0.15, wspace=0.15)

    # Compute weighted correlation matrix
    X = df[PARAMS].to_numpy(float)
    mean_w = np.average(X, weights=w, axis=0)
    Xc = X - mean_w
    cov_w = np.dot(w * Xc.T, Xc)
    diag_std = np.sqrt(np.diag(cov_w))
    corr_mat = cov_w / np.outer(diag_std, diag_std)

    for i, row_p in enumerate(PARAMS):
        for j, col_p in enumerate(PARAMS):
            ax = axes[i, j]

            # Diagonal: 1D marginal posterior KDE vs Prior
            if i == j:
                vals = df[row_p].to_numpy(float)
                p_lo, p_hi = PRIORS[row_p]
                grid_x = np.linspace(p_lo, p_hi, 200)

                kde = stats.gaussian_kde(vals, weights=w, bw_method="scott")
                dens = kde(grid_x)

                ax.plot(grid_x, dens, color="#1D4ED8", lw=1.6, label="Posterior KDE")
                ax.fill_between(grid_x, 0, dens, color="#93C5FD", alpha=0.35)

                prior_dens = 1.0 / (p_hi - p_lo)
                ax.axhline(prior_dens, color="#9CA3AF", linestyle="--", lw=1.2, label="Uniform prior")

                med, q025, q975 = weighted_stats(vals, w)
                ax.axvline(med, color="#DC2626", linestyle="-", lw=1.0)

                ax.set_xlim(p_lo, p_hi)
                ax.set_ylim(bottom=0)
                ax.set_xticks(TICKS[row_p])
                ax.tick_params(labelsize=8, length=2.5)

                if i == 0:
                    ax.set_title(PARAM_LABELS[row_p], fontsize=11, fontweight="bold", pad=8)
                else:
                    ax.set_title("")

                if j == 0:
                    ax.set_ylabel(f"{PARAM_LABELS[row_p]}\nDensity", fontsize=9.5)
                else:
                    ax.set_ylabel("")

            # Upper triangle: Scatter points + linear regression trendline + r annotation
            elif i < j:
                x_vals = df[col_p].to_numpy(float)
                y_vals = df[row_p].to_numpy(float)
                r_val = corr_mat[i, j]

                # Scatter points sized and shaded by ABC particle weight
                norm_w = (w - w.min()) / (w.max() - w.min() + 1e-12)
                point_sizes = 15 + 40 * norm_w
                ax.scatter(x_vals, y_vals, s=point_sizes, color="#2563EB", alpha=0.55, edgecolors="#1E3A8A", linewidths=0.4)

                # Linear trend line
                poly = np.polyfit(x_vals, y_vals, 1, w=np.sqrt(w))
                x_trend = np.linspace(PRIORS[col_p][0], PRIORS[col_p][1], 50)
                y_trend = poly[0] * x_trend + poly[1]
                ax.plot(x_trend, y_trend, color="#DC2626", linestyle="--", lw=1.2)

                # Correlation label
                ax.text(
                    0.07, 0.88,
                    f"$r = {r_val:+.2f}$",
                    transform=ax.transAxes,
                    fontsize=9.5,
                    fontweight="bold",
                    color="#B91C1C" if abs(r_val) >= 0.4 else "#374151",
                    bbox=dict(boxstyle="round,pad=0.2", facecolor="#F3F4F6", edgecolor="none", alpha=0.8),
                )

                ax.set_xlim(PRIORS[col_p])
                ax.set_ylim(PRIORS[row_p])
                ax.set_xticks(TICKS[col_p])
                ax.set_yticks(TICKS[row_p])
                ax.tick_params(labelsize=8, length=2.5)

                if i == 0:
                    ax.set_title(PARAM_LABELS[col_p], fontsize=11, fontweight="bold", pad=8)

            # Lower triangle: 2D Bivariate KDE contour density
            else:
                x_vals = df[col_p].to_numpy(float)
                y_vals = df[row_p].to_numpy(float)
                r_val = corr_mat[i, j]

                # Bivariate KDE
                p_lo_x, p_hi_x = PRIORS[col_p]
                p_lo_y, p_hi_y = PRIORS[row_p]
                gx = np.linspace(p_lo_x, p_hi_x, 80)
                gy = np.linspace(p_lo_y, p_hi_y, 80)
                GX, GY = np.meshgrid(gx, gy)
                positions = np.vstack([GX.ravel(), GY.ravel()])
                values = np.vstack([x_vals, y_vals])

                kde2d = stats.gaussian_kde(values, weights=w, bw_method="scott")
                GZ = np.reshape(kde2d(positions).T, GX.shape)

                # Plot contour levels
                ax.contourf(GX, GY, GZ, levels=7, cmap="YlOrRd", alpha=0.75)
                ax.contour(GX, GY, GZ, levels=7, colors="#7C2D12", linewidths=0.5, alpha=0.7)

                # Overlay individual particles
                ax.scatter(x_vals, y_vals, s=8, color="#1F2937", alpha=0.35, edgecolors="none")

                # Correlation label
                ax.text(
                    0.07, 0.88,
                    f"$r = {r_val:+.2f}$",
                    transform=ax.transAxes,
                    fontsize=9.5,
                    fontweight="bold",
                    color="#B91C1C" if abs(r_val) >= 0.4 else "#374151",
                    bbox=dict(boxstyle="round,pad=0.2", facecolor="#FFFFFF", edgecolor="#E5E7EB", alpha=0.85),
                )

                ax.set_xlim(p_lo_x, p_hi_x)
                ax.set_ylim(p_lo_y, p_hi_y)
                ax.set_xticks(TICKS[col_p])
                ax.set_yticks(TICKS[row_p])
                ax.tick_params(labelsize=8, length=2.5)

            # Labels only on bottom row and left column
            if i == n_params - 1:
                ax.set_xlabel(PARAM_LABELS[col_p], fontsize=10.5, fontweight="bold")
            else:
                ax.set_xticklabels([])

            if j == 0:
                ax.set_ylabel(PARAM_LABELS[row_p], fontsize=10.5, fontweight="bold")
            else:
                ax.set_yticklabels([])

    fig.suptitle("ABC-SMC Calibrated Parameter Pairwise Scatter Matrix (Homeostasis Model)", fontsize=12, fontweight="bold", y=0.99)

    out_dir.mkdir(parents=True, exist_ok=True)
    pdf_path = out_dir / "Figure4_pairwise_scatters.pdf"
    png_path = out_dir / "Figure4_pairwise_scatters.png"

    plt.savefig(pdf_path, format="pdf", bbox_inches="tight")
    plt.savefig(png_path, format="png", dpi=300, bbox_inches="tight")
    plt.close()

    print(f"Saved: {pdf_path}")
    print(f"Saved: {png_path}")
    return pdf_path, png_path


def main():
    parser = argparse.ArgumentParser(description="Generate publication pairwise scatter matrix for Figure 4 ABC parameters.")
    parser.add_argument(
        "--data-file",
        type=Path,
        default=Path("experiments/exp08_cs_abc_smc_homeostasis/data/generation_14.csv"),
        help="Path to posterior particles CSV file",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("experiments/exp08_cs_abc_smc_homeostasis"),
        help="Output directory",
    )
    args = parser.parse_args()

    pdf_path, png_path = render_pairwise_scatters(args.data_file, args.output_dir)

    # Mirror to results/paper_figures
    paper_dir = Path("results/paper_figures")
    paper_dir.mkdir(parents=True, exist_ok=True)
    import shutil
    shutil.copy2(pdf_path, paper_dir / "Figure4_pairwise_scatters.pdf")
    shutil.copy2(png_path, paper_dir / "Figure4_pairwise_scatters.png")
    print(f"Mirrored to: {paper_dir / 'Figure4_pairwise_scatters.png'}")


if __name__ == "__main__":
    main()
