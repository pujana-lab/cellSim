#!/usr/bin/env python3
"""
generate_sfig4_1000cells.py — Supplementary Figures S4 & S4b for 1000-cell ABC-SMC calibration.

Adapted from experiments/exp03_cs_abc_smc_calibration/plot_sfig4.py:
  - Figure S4 (Variant 4a): Corner matrix (4x4) for the 4 canonical parameters
    (xi, low_delta, high_delta, neoplastic_div_rate/nu).
  - Figure S4b (Variant 4b): 2x3 Grid of 1D posterior distributions (SIS weighted KDE)
    and 95% CrIs for all 6 pairwise parameter ratios.
"""

from __future__ import annotations

import argparse
import shutil
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Tuple

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


@dataclass
class VariantConfig:
    title: str
    file_stem: str
    params: List[str]
    labels: Dict[str, str]
    axis_limits: Dict[str, Tuple[float, float]]
    ticks: Dict[str, List[float]]


def weighted_stats(x: np.ndarray, w: np.ndarray) -> Tuple[float, float, float]:
    """Calculate weighted median, 2.5%, and 97.5% percentiles."""
    idx = np.argsort(x)
    sx = x[idx]
    sw = w[idx]
    cw = np.cumsum(sw) / np.sum(sw)
    q025 = float(np.interp(0.025, cw, sx))
    q50 = float(np.interp(0.50, cw, sx))
    q975 = float(np.interp(0.975, cw, sx))
    return q50, q025, q975


def load_particles(data_csv: Path, use_xi: bool = True) -> Tuple[pd.DataFrame, np.ndarray]:
    """Load posterior particles and compute all 6 pairwise parameter ratios."""
    if not data_csv.exists():
        raise FileNotFoundError(f"Missing {data_csv}")

    df = pd.read_csv(data_csv)
    w = df["weight"].to_numpy(float)
    w = w / np.sum(w)

    dl = df["low_delta"]
    dh = df["high_delta"]
    nu = df["neoplastic_div_rate"]

    if use_xi:
        mut_var = df["xi"]
        df["ratio_dh_dl"] = dh / dl
        df["ratio_dl_mut"] = dl / mut_var
        df["ratio_nu_mut"] = nu / mut_var
        df["ratio_nu_dl"] = nu / dl
        df["ratio_dh_mut"] = dh / mut_var
        df["ratio_nu_dh"] = nu / dh
    else:
        mut_var = df["beta"]
        df["ratio_dh_dl"] = dh / dl
        df["ratio_dl_mut"] = dl / mut_var
        df["ratio_nu_mut"] = nu / mut_var
        df["ratio_nu_dl"] = nu / dl
        df["ratio_dh_mut"] = dh / mut_var
        df["ratio_nu_dh"] = nu / dh

    return df, w


def render_corner_matrix_s4(df: pd.DataFrame, w: np.ndarray, cfg: VariantConfig, out_dir: Path) -> Tuple[Path, Path]:
    """Renders the 4x4 posterior corner matrix (Figure S4)."""
    n = len(cfg.params)
    fig, axes = plt.subplots(n, n, figsize=(8.8, 8.8), dpi=300)

    fig.text(0.06, 0.98, "Supplementary Figure S4", fontsize=15, fontweight="bold", va="top")
    fig.text(0.06, 0.95, "Joint and marginal posterior parameter distributions (1000 cells neutral dynamics, F = 1.0)",
             fontsize=10.0, color="#374151", va="top")

    for i, p_row in enumerate(cfg.params):
        for j, p_col in enumerate(cfg.params):
            ax = axes[i, j]

            # Diagonal: 1D KDE profile
            if i == j:
                vals = df[p_row].to_numpy(float)
                xlim = cfg.axis_limits[p_row]
                grid_x = np.linspace(xlim[0], xlim[1], 300)

                kde = stats.gaussian_kde(vals, weights=w)
                dens = kde(grid_x)
                dens_max = np.max(dens) if np.max(dens) > 0 else 1.0
                dens_norm = dens / dens_max

                ax.fill_between(grid_x, 0, dens_norm, color="#FED7AA", alpha=0.85)
                ax.plot(grid_x, dens_norm, color="#7C2D12", linewidth=1.1)

                med, q025, q975 = weighted_stats(vals, w)
                ax.axvline(med, color="#C2410C", linewidth=1.2)
                ax.axvline(q025, color="#C2410C", linestyle="--", linewidth=1.0)
                ax.axvline(q975, color="#C2410C", linestyle="--", linewidth=1.0)

                ax.set_xlim(xlim)
                ax.set_ylim(0.0, 1.05)
                if p_row in cfg.ticks:
                    ax.set_xticks(cfg.ticks[p_row])
                ax.set_yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
                ax.set_title(cfg.labels[p_row], fontsize=12, pad=4)
                ax.tick_params(labelsize=8, direction="out", length=2.5)

                if i == 0 and j == 0:
                    ax.set_ylabel(cfg.labels[p_row], fontsize=12)
                else:
                    if j > 0:
                        ax.set_yticklabels([])

                if i == n - 1:
                    ax.set_xlabel(cfg.labels[p_col], fontsize=12)
                else:
                    ax.set_xticklabels([])

            # Lower Triangle: 2D Bivariate Density Contours
            elif i > j:
                x = df[p_col].to_numpy(float)
                y = df[p_row].to_numpy(float)
                xlim = cfg.axis_limits[p_col]
                ylim = cfg.axis_limits[p_row]

                kde2d = stats.gaussian_kde(np.vstack([x, y]), weights=w)
                gx = np.linspace(xlim[0], xlim[1], 120)
                gy = np.linspace(ylim[0], ylim[1], 120)
                X, Y = np.meshgrid(gx, gy)
                Z = kde2d(np.vstack([X.ravel(), Y.ravel()])).reshape(X.shape)
                z_max = np.max(Z) if np.max(Z) > 0 else 1.0
                Z_norm = Z / z_max

                levels = [0.08, 0.22, 0.42, 0.65, 0.88, 1.0]
                ax.contourf(X, Y, Z_norm, levels=levels, cmap="YlOrBr", alpha=0.9)
                ax.contour(X, Y, Z_norm, levels=levels, colors="#78350F", linewidths=0.7, alpha=0.9)
                ax.scatter(x, y, s=6, color="#451A03", alpha=0.35, edgecolors="none")

                ax.set_xlim(xlim)
                ax.set_ylim(ylim)
                if p_col in cfg.ticks:
                    ax.set_xticks(cfg.ticks[p_col])
                if p_row in cfg.ticks:
                    ax.set_yticks(cfg.ticks[p_row])
                ax.tick_params(labelsize=8, direction="out", length=2.5)

                if j == 0:
                    ax.set_ylabel(cfg.labels[p_row], fontsize=12)
                else:
                    ax.set_yticklabels([])

                if i == n - 1:
                    ax.set_xlabel(cfg.labels[p_col], fontsize=12)
                else:
                    ax.set_xticklabels([])

            # Upper Triangle: Color Gradient Legend Box
            else:
                ax.axis("off")
                if i == 0 and j == 1:
                    rect = plt.Rectangle((0.05, 0.15), 0.90, 0.75, transform=ax.transAxes,
                                         facecolor="white", edgecolor="#111827", linewidth=0.8)
                    ax.add_patch(rect)
                    ax.text(0.10, 0.78, "Normalized posterior density", transform=ax.transAxes,
                            fontsize=9, va="top", color="#111827")

                    cax = ax.inset_axes([0.15, 0.25, 0.10, 0.42])
                    norm = mpl.colors.Normalize(vmin=0.0, vmax=1.0)
                    cb = fig.colorbar(
                        mpl.cm.ScalarMappable(norm=norm, cmap="YlOrBr"),
                        cax=cax, ticks=[0.0, 0.5, 1.0]
                    )
                    cb.ax.tick_params(labelsize=7.5, length=2.0)
                    cb.outline.set_linewidth(0.6)

    fig.subplots_adjust(left=0.10, right=0.95, bottom=0.10, top=0.92, hspace=0.18, wspace=0.18)

    out_dir.mkdir(parents=True, exist_ok=True)
    pdf_path = out_dir / f"{cfg.file_stem}.pdf"
    png_path = out_dir / f"{cfg.file_stem}.png"
    fig.savefig(pdf_path, bbox_inches="tight")
    fig.savefig(png_path, dpi=300, bbox_inches="tight")
    plt.close(fig)

    print(f"Saved: {pdf_path}")
    print(f"Saved: {png_path}")
    return pdf_path, png_path


def render_sfig4b_ratios(df: pd.DataFrame, w: np.ndarray, out_dir: Path, use_xi: bool = True) -> Tuple[Path, Path]:
    """Renders Figure S4b: 2x3 grid displaying 1D posterior distributions of pairwise ratios."""
    mut_sym = r"\xi" if use_xi else r"\beta"
    mut_desc = "Mutation Ratio" if use_xi else "Mutation Rate"

    ratio_specs = [
        ("ratio_dh_dl", r"$\delta_{\rm high} / \delta_{\rm low}$", "Accelerated vs Latent Accumulation", "a"),
        ("ratio_dl_mut", rf"$\delta_{{\rm low}} / {mut_sym}$", f"Latent Accumulation vs {mut_desc}", "b"),
        ("ratio_nu_mut", rf"$\nu / {mut_sym}$", f"Tumoral Proliferation vs {mut_desc}", "c"),
        ("ratio_nu_dl", r"$\nu / \delta_{\rm low}$", "Tumoral Proliferation vs Latent Accumulation", "d"),
        ("ratio_dh_mut", rf"$\delta_{{\rm high}} / {mut_sym}$", f"Accelerated Accumulation vs {mut_desc}", "e"),
        ("ratio_nu_dh", r"$\nu / \delta_{\rm high}$", "Tumoral Proliferation vs Accelerated Accumulation", "f"),
    ]

    fig, axes = plt.subplots(2, 3, figsize=(12.0, 7.8), dpi=300)

    fig.text(0.06, 0.985, "Supplementary Figure S4b", fontsize=14, fontweight="bold", va="top")
    fig.text(0.06, 0.948,
             r"Posterior probability densities of all pairwise parameter ratios ($R = \theta_i / \theta_j$) under 1000 cells neutral dynamics",
             fontsize=10.0, color="#374151", va="top")

    for idx, (col_name, ratio_lbl, desc, letter) in enumerate(ratio_specs):
        row, col = divmod(idx, 3)
        ax = axes[row, col]

        vals = df[col_name].to_numpy(float)
        med, q025, q975 = weighted_stats(vals, w)

        val_min, val_max = np.min(vals), np.max(vals)
        val_span = val_max - val_min
        x_min = max(0.0, val_min - 0.15 * val_span)
        x_max = val_max + 0.15 * val_span

        grid_x = np.linspace(x_min, x_max, 300)
        kde = stats.gaussian_kde(vals, weights=w)
        dens = kde(grid_x)

        cri_mask = (grid_x >= q025) & (grid_x <= q975)

        ax.fill_between(grid_x, 0, dens, color="#FED7AA", alpha=0.50, label="Posterior density")
        ax.fill_between(grid_x[cri_mask], 0, dens[cri_mask], color="#FDBA74", alpha=0.85, label="95% Credible interval")
        ax.plot(grid_x, dens, color="#7C2D12", linewidth=1.4)

        ax.scatter(vals, np.zeros_like(vals), color="#451A03", alpha=0.35, s=12, marker="|", zorder=3)

        # Reference lines
        label_med = f"Median ({med:.3f})" if med < 0.1 else f"Median ({med:.2f})"
        ax.axvline(med, color="#C2410C", linewidth=1.6, label=label_med)
        ax.axvline(q025, color="#9A3412", linestyle="--", linewidth=1.0)
        ax.axvline(q975, color="#9A3412", linestyle="--", linewidth=1.0)

        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.set_xlim(x_min, x_max)
        ax.set_ylim(bottom=0.0, top=np.max(dens) * 1.25)
        ax.set_xlabel(f"Ratio {ratio_lbl}", fontsize=10, labelpad=4)
        if col == 0:
            ax.set_ylabel("Posterior density", fontsize=10, labelpad=4)
        else:
            ax.set_ylabel("")

        ax.set_title(f"({letter}) {ratio_lbl}\n{desc}", fontsize=9.5, fontweight="bold", pad=6)
        ax.tick_params(labelsize=8.5, direction="out", length=2.5)

        annot_text = (
            f"Median: {med:.3f}\n95% CrI: [{q025:.3f}, {q975:.3f}]" if med < 0.1
            else f"Median: {med:.2f}\n95% CrI: [{q025:.2f}, {q975:.2f}]"
        )
        ax.text(
            0.96, 0.94, annot_text, transform=ax.transAxes,
            fontsize=8.5, ha="right", va="top",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor="#D1D5DB", alpha=0.92)
        )

        if idx == 0:
            ax.legend(loc="upper left", fontsize=7.0, frameon=False)

    fig.subplots_adjust(left=0.07, right=0.96, bottom=0.08, top=0.84, hspace=0.45, wspace=0.22)

    out_dir.mkdir(parents=True, exist_ok=True)
    pdf_path = out_dir / "FigureS4b.pdf"
    png_path = out_dir / "FigureS4b.png"
    fig.savefig(pdf_path, bbox_inches="tight")
    fig.savefig(png_path, dpi=300, bbox_inches="tight")
    plt.close(fig)

    print(f"Saved: {pdf_path}")
    print(f"Saved: {png_path}")
    return pdf_path, png_path


def main():
    parser = argparse.ArgumentParser(description="Generate Figure S4 and S4b for 1000 cells calibration.")
    parser.add_argument(
        "--data-file",
        type=Path,
        default=Path("experiments/exp_101_cs_xi_contour/runs/run_abc_xi_wide_20261120_171832/particles_gen_10.csv"),
        help="Path to posterior particles CSV",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("experiments/exp_101_cs_xi_contour/results_1000cells_final"),
        help="Output directory",
    )
    parser.add_argument(
        "--param",
        choices=["xi", "beta"],
        default="xi",
        help="Whether to use xi or beta as the primary mutation parameter (default: xi)",
    )
    args = parser.parse_args()

    use_xi = (args.param == "xi")
    df, w = load_particles(args.data_file, use_xi=use_xi)

    if use_xi:
        cfg = VariantConfig(
            title="Figure S4",
            file_stem="FigureS4",
            params=["xi", "low_delta", "high_delta", "neoplastic_div_rate"],
            labels={
                "xi": r"$\xi$",
                "low_delta": r"$\delta_{\rm low}$",
                "high_delta": r"$\delta_{\rm high}$",
                "neoplastic_div_rate": r"$\nu$",
            },
            axis_limits={
                "xi": (10.0, 52.0),
                "low_delta": (0.09, 0.27),
                "high_delta": (0.14, 0.54),
                "neoplastic_div_rate": (0.18, 0.56),
            },
            ticks={
                "xi": [20.0, 30.0, 40.0, 50.0],
                "low_delta": [0.10, 0.15, 0.20, 0.25],
                "high_delta": [0.20, 0.35, 0.50],
                "neoplastic_div_rate": [0.20, 0.35, 0.50],
            },
        )
    else:
        cfg = VariantConfig(
            title="Figure S4",
            file_stem="FigureS4",
            params=["beta", "low_delta", "high_delta", "neoplastic_div_rate"],
            labels={
                "beta": r"$\beta$",
                "low_delta": r"$\delta_{\rm low}$",
                "high_delta": r"$\delta_{\rm high}$",
                "neoplastic_div_rate": r"$\nu$",
            },
            axis_limits={
                "beta": (0.005, 0.075),
                "low_delta": (0.09, 0.27),
                "high_delta": (0.14, 0.54),
                "neoplastic_div_rate": (0.18, 0.56),
            },
            ticks={
                "beta": [0.01, 0.03, 0.05, 0.07],
                "low_delta": [0.10, 0.15, 0.20, 0.25],
                "high_delta": [0.20, 0.35, 0.50],
                "neoplastic_div_rate": [0.20, 0.35, 0.50],
            },
        )

    s4_pdf, s4_png = render_corner_matrix_s4(df, w, cfg, args.output_dir)
    s4b_pdf, s4b_png = render_sfig4b_ratios(df, w, args.output_dir, use_xi=use_xi)

    print(f"\n[DONE] Generadas Figure S4 y Figure S4b (con {args.param}) en: {args.output_dir}")


if __name__ == "__main__":
    main()
