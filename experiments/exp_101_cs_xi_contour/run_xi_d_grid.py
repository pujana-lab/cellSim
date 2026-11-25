#!/usr/bin/env python3
"""
run_xi_d_grid.py — Evaluates penetrance on a 2D logarithmic grid over:
  d in [10^-3, 10^-2]
  xi in [1, 100]
With exact analytical neutrality / homeostasis F = 1.0 (lambda_B = 1.0):
  c = F^(1/T) / (1 + d)
  mu = 2 * (1 - c) / ((1 + xi) + sqrt((1 + xi)^2 - 4 * xi * (1 - c)))
  beta = xi * mu

Calculates cumulative penetrance at ages 30, 40, 50, 60, 70, 80 years
and compares against Kuchenbaecker 2017 to produce the 2D contour map requested by Víctor.
"""

from __future__ import annotations

import os
import re
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

plt.switch_backend("Agg")

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parents[1]
DATA_DIR = SCRIPT_DIR / "data"
FIG_DIR = SCRIPT_DIR / "figures"
BOOTSTRAP_BIN = PROJECT_ROOT / "build" / "run_bootstrapping"


def mu_beta(d: float, xi: float, F: float = 1.0, T: int = 80) -> tuple[float, float]:
    c = (F ** (1.0 / T)) / (1.0 + d)
    term = (1.0 + xi) ** 2 - 4.0 * xi * (1.0 - c)
    if term < 0:
        raise ValueError(f"Negative discriminant for d={d}, xi={xi}, F={F}")
    mu = 2.0 * (1.0 - c) / ((1.0 + xi) + np.sqrt(term))
    beta = xi * mu
    assert np.allclose(((1.0 - beta) * (1.0 - mu) * (1.0 + d)) ** T, F, atol=1e-5)
    return float(mu), float(beta)


def run_single_point(args: tuple) -> dict:
    d_val, xi_val, F_val, n_runs, n_cells, delta_low, delta_high, seed_offset = args
    mu_val, beta_val = mu_beta(d_val, xi_val, F_val)

    cmd = [
        str(BOOTSTRAP_BIN),
        "--n-runs", str(n_runs),
        "--weighted-sse",
        "--big-bang",
        "--n-cells", str(n_cells),
        "--max-healthy-cells", str(n_cells),
        "--division-rate", f"{d_val:.8f}",
        "--tp53-rate", f"{mu_val:.8f}",
        "--brca1-rate", f"{beta_val:.8f}",
        "--low-delta", str(delta_low),
        "--high-delta", str(delta_high),
        "--neoplastic-div-rate", "0.250",
        "--output", "/dev/null",
        "--seed-offset", str(seed_offset),
    ]

    env = dict(os.environ)
    env["OMP_NUM_THREADS"] = "1"
    res = subprocess.run(cmd, capture_output=True, text=True, env=env)

    sse_match = re.search(r"weighted SSE\s*=\s*([\d.eE+\-]+)", res.stdout)
    sse_val = float(sse_match.group(1)) if sse_match else np.nan

    risks = {}
    for m in re.finditer(r"^\s*(\d{2})\s+([\d.]+)\s+", res.stdout, re.MULTILINE):
        risks[int(m.group(1))] = float(m.group(2))

    return {
        "d": d_val,
        "xi": xi_val,
        "mu": mu_val,
        "beta": beta_val,
        "F": F_val,
        "weighted_sse": sse_val,
        "risk_30": risks.get(30, np.nan),
        "risk_40": risks.get(40, np.nan),
        "risk_50": risks.get(50, np.nan),
        "risk_60": risks.get(60, np.nan),
        "risk_70": risks.get(70, np.nan),
        "risk_80": risks.get(80, np.nan),
    }


def main():
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    # 1. Configuración del Grid logarítmico (8 x 8 = 64 puntos)
    n_d = 8
    n_xi = 8
    d_vals = np.logspace(-3, -2, n_d)  # 10^-3 a 10^-2
    xi_vals = np.logspace(0, 2, n_xi)   # 1 a 100

    F_target = 1.0
    n_runs = 250
    n_cells = 500
    delta_low = 0.102
    delta_high = 0.200

    tasks = []
    seed_base = 10000
    idx = 0
    for d in d_vals:
        for xi in xi_vals:
            tasks.append((d, xi, F_target, n_runs, n_cells, delta_low, delta_high, seed_base + idx * 50))
            idx += 1

    print("=" * 70)
    print("EJECUTANDO GRID BIDIMENSIONAL: d in [10^-3, 10^-2] x XI in [1, 100]")
    print(f"Dimensiones: {n_d} x {n_xi} = {len(tasks)} puntos de simulación")
    print(f"Condición: F = {F_target} (Homeostasis neutra lambda_B = 1.0)")
    print(f"Células: {n_cells} | Réplicas por punto: {n_runs}")
    print("=" * 70)

    t0 = time.time()
    with ProcessPoolExecutor(max_workers=10) as executor:
        results = list(executor.map(run_single_point, tasks))
    t1 = time.time()
    print(f"Simulaciones finalizadas en {t1 - t0:.1f} segundos.")

    df = pd.DataFrame(results)
    csv_out = DATA_DIR / "grid_xi_d_penetrance.csv"
    df.to_csv(csv_out, index=False)
    print(f"Datos guardados en: {csv_out}")

    # 2. Generar Figuras de Contorno estilo Víctor
    plot_contour_figure(df, d_vals, xi_vals)


def plot_contour_figure(df: pd.DataFrame, d_vals: np.ndarray, xi_vals: np.ndarray):
    # Pivot matrices
    grid_r80 = df.pivot(index="xi", columns="d", values="risk_80").values
    grid_r50 = df.pivot(index="xi", columns="d", values="risk_50").values
    grid_sse = df.pivot(index="xi", columns="d", values="weighted_sse").values

    fig, axes = plt.subplots(1, 3, figsize=(18, 5.5), dpi=300)

    # Panel 1: Penetrancia a los 80 años con contornos
    ax1 = axes[0]
    im1 = ax1.pcolormesh(d_vals, xi_vals, grid_r80, shading="gouraud", cmap="Blues", vmin=0, vmax=100)
    c1 = ax1.contour(d_vals, xi_vals, grid_r80, levels=[40, 60, 70, 80, 90, 95], colors="white", linewidths=1.2)
    ax1.clabel(c1, inline=True, fontsize=8, fmt="%1.0f%%")
    ax1.set_xscale("log")
    ax1.set_yscale("log")
    ax1.set_title("a   Penetrance at age 80 (%)", fontsize=12, fontweight="bold")
    ax1.set_xlabel("Division rate $d$ (log scale)", fontsize=11)
    ax1.set_ylabel(r"Ratio $\xi = \beta / \mu$ (log scale)", fontsize=11)
    fig.colorbar(im1, ax=ax1, label="Penetrance at 80y (%)")

    # Panel 2: Superposición de contornos Edad 80 (solid) vs Edad 50 (dashed) como en Panel D de Víctor
    ax2 = axes[1]
    im2 = ax2.pcolormesh(d_vals, xi_vals, grid_r80, shading="gouraud", cmap="Blues", vmin=0, vmax=100, alpha=0.85)
    cs80 = ax2.contour(d_vals, xi_vals, grid_r80, levels=[50, 70, 85, 95], colors="black", linestyles="solid", linewidths=1.5)
    cs50 = ax2.contour(d_vals, xi_vals, grid_r50, levels=[20, 35, 46, 60], colors="crimson", linestyles="dashed", linewidths=1.5)
    ax2.clabel(cs80, inline=True, fontsize=8, fmt="80y: %1.0f%%")
    ax2.clabel(cs50, inline=True, fontsize=8, fmt="50y: %1.0f%%")
    ax2.set_xscale("log")
    ax2.set_yscale("log")
    ax2.set_title("b   Contour Shift: Age 80 (solid) vs Age 50 (dashed)", fontsize=12, fontweight="bold")
    ax2.set_xlabel("Division rate $d$ (log scale)", fontsize=11)
    ax2.set_ylabel(r"Ratio $\xi = \beta / \mu$ (log scale)", fontsize=11)
    # Leyenda manual
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], color="black", lw=1.5, linestyle="solid", label="Age 80 contours"),
        Line2D([0], [0], color="crimson", lw=1.5, linestyle="dashed", label="Age 50 contours (46% clinical)"),
    ]
    ax2.legend(handles=legend_elements, loc="lower right", framealpha=0.9, fontsize=9)

    # Panel 3: Weighted SSE frente a Kuchenbaecker 2017
    ax3 = axes[2]
    im3 = ax3.pcolormesh(d_vals, xi_vals, grid_sse, shading="gouraud", cmap="viridis_r")
    cs3 = ax3.contour(d_vals, xi_vals, grid_sse, levels=8, colors="white", linewidths=1.0)
    ax3.clabel(cs3, inline=True, fontsize=8, fmt="%1.0f")
    ax3.set_xscale("log")
    ax3.set_yscale("log")
    ax3.set_title("c   Weighted SSE vs Kuchenbaecker 2017", fontsize=12, fontweight="bold")
    ax3.set_xlabel("Division rate $d$ (log scale)", fontsize=11)
    ax3.set_ylabel(r"Ratio $\xi = \beta / \mu$ (log scale)", fontsize=11)
    fig.colorbar(im3, ax=ax3, label="Weighted SSE (lower is better)")

    plt.tight_layout()
    out_png = FIG_DIR / "Figure_grid_xi_d_contour.png"
    out_pdf = FIG_DIR / "Figure_grid_xi_d_contour.pdf"
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.savefig(out_pdf, bbox_inches="tight")
    plt.close()

    print(f"\n[DONE] Gráficas de contorno generadas:")
    print(f"  -> PNG: {out_png}")
    print(f"  -> PDF: {out_pdf}")


if __name__ == "__main__":
    main()
