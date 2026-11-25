#!/usr/bin/env python3
"""
run_xi_d_compare_F.py — Runs fast 2D grids for F = 0.75 and F = 1.25.
Uses N=100 Monte Carlo runs (faster execution, ~1-2 min total).
"""

from __future__ import annotations

import os
import re
import subprocess
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


def mu_beta(d: float, xi: float, F: float, T: int = 80) -> tuple[float, float] | tuple[None, None]:
    c = (F ** (1.0 / T)) / (1.0 + d)
    if c >= 1.0:
        # Imposible biológicamente: las mutaciones beta, mu tendrían que ser negativas
        return None, None
    term = (1.0 + xi) ** 2 - 4.0 * xi * (1.0 - c)
    if term < 0:
        return None, None
    mu = 2.0 * (1.0 - c) / ((1.0 + xi) + np.sqrt(term))
    beta = xi * mu
    if mu < 0 or beta < 0:
        return None, None
    return float(mu), float(beta)


def run_single_point(args: tuple) -> dict:
    d_val, xi_val, F_val, n_runs, n_cells, delta_low, delta_high, seed_offset = args
    mu_val, beta_val = mu_beta(d_val, xi_val, F_val)

    if mu_val is None or beta_val is None:
        return {
            "d": d_val, "xi": xi_val, "mu": np.nan, "beta": np.nan, "F": F_val,
            "weighted_sse": np.nan, "risk_50": np.nan, "risk_80": np.nan, "valid": False
        }

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
        "risk_50": risks.get(50, np.nan),
        "risk_80": risks.get(80, np.nan),
        "valid": True,
    }


def main():
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    # Grilla rápida de 6x6 (36 puntos por valor de F)
    n_d = 6
    n_xi = 6
    d_vals = np.logspace(-3, -2, n_d)  # 10^-3 a 10^-2
    xi_vals = np.logspace(0, 2, n_xi)   # 1 a 100

    n_runs = 100
    n_cells = 500
    delta_low = 0.102
    delta_high = 0.200

    for F_target in [0.75, 1.25]:
        print("=" * 65)
        print(f"LANZANDO GRID RÁPIDO PARA F = {F_target} (N={n_runs} runs, 6x6 grid)")
        print("=" * 65)

        tasks = []
        seed_base = int(F_target * 100000)
        idx = 0
        for d in d_vals:
            for xi in xi_vals:
                tasks.append((d, xi, F_target, n_runs, n_cells, delta_low, delta_high, seed_base + idx * 50))
                idx += 1

        t0 = time.time()
        with ProcessPoolExecutor(max_workers=10) as executor:
            results = list(executor.map(run_single_point, tasks))
        t1 = time.time()
        print(f"F={F_target} completado en {t1 - t0:.1f} segundos.")

        df = pd.DataFrame(results)
        csv_path = DATA_DIR / f"grid_xi_d_F_{F_target:.2f}.csv"
        df.to_csv(csv_path, index=False)

    # Generar figura comparativa de F=0.75 vs F=1.00 vs F=1.25
    plot_comparative_F(d_vals, xi_vals)


def plot_comparative_F(d_vals, xi_vals):
    df_075 = pd.read_csv(DATA_DIR / "grid_xi_d_F_0.75.csv")
    df_125 = pd.read_csv(DATA_DIR / "grid_xi_d_F_1.25.csv")

    fig, axes = plt.subplots(2, 2, figsize=(14, 11), dpi=300)

    # F=0.75: Riesgo 80 años y Contornos
    ax1 = axes[0, 0]
    grid_075_80 = df_075.pivot(index="xi", columns="d", values="risk_80").values
    im1 = ax1.pcolormesh(d_vals, xi_vals, grid_075_80, shading="gouraud", cmap="Blues", vmin=0, vmax=100)
    c1 = ax1.contour(d_vals, xi_vals, grid_075_80, levels=[30, 50, 70, 90], colors="white", linewidths=1.2)
    ax1.clabel(c1, inline=True, fontsize=8, fmt="%1.0f%%")
    ax1.set_xscale("log")
    ax1.set_yscale("log")
    ax1.set_title("a   F = 0.75 (Contracción tisular 25%): Penetrancia a 80 años", fontsize=11, fontweight="bold")
    ax1.set_xlabel("d (log scale)")
    ax1.set_ylabel(r"$\xi = \beta/\mu$ (log scale)")
    fig.colorbar(im1, ax=ax1)

    # F=0.75: Shift 50 vs 80
    ax2 = axes[0, 1]
    grid_075_50 = df_075.pivot(index="xi", columns="d", values="risk_50").values
    ax2.pcolormesh(d_vals, xi_vals, grid_075_80, shading="gouraud", cmap="Blues", vmin=0, vmax=100, alpha=0.8)
    cs80_1 = ax2.contour(d_vals, xi_vals, grid_075_80, levels=[50, 70, 90], colors="black", linestyles="solid")
    cs50_1 = ax2.contour(d_vals, xi_vals, grid_075_50, levels=[20, 46, 60], colors="crimson", linestyles="dashed")
    ax2.clabel(cs80_1, inline=True, fontsize=8, fmt="80y: %1.0f%%")
    ax2.clabel(cs50_1, inline=True, fontsize=8, fmt="50y: %1.0f%%")
    ax2.set_xscale("log")
    ax2.set_yscale("log")
    ax2.set_title("b   F = 0.75: Desplazamiento 80a (sólido) vs 50a (discontinuo)", fontsize=11, fontweight="bold")
    ax2.set_xlabel("d (log scale)")
    ax2.set_ylabel(r"$\xi = \beta/\mu$ (log scale)")

    # F=1.25: Riesgo 80 años
    ax3 = axes[1, 0]
    grid_125_80 = df_125.pivot(index="xi", columns="d", values="risk_80").values
    im3 = ax3.pcolormesh(d_vals, xi_vals, grid_125_80, shading="gouraud", cmap="Blues", vmin=0, vmax=100)
    # Mascarar los puntos matemáticamente no válidos (donde d no da abasto para F=1.25)
    valid_mask = df_125.pivot(index="xi", columns="d", values="valid").values
    for i in range(len(xi_vals)):
        for j in range(len(d_vals)):
            if not valid_mask[i, j]:
                ax3.text(d_vals[j], xi_vals[i], "N/A", ha="center", va="center", color="gray", fontsize=7)
    c3 = ax3.contour(d_vals, xi_vals, grid_125_80, levels=[30, 50, 70, 90], colors="white", linewidths=1.2)
    ax3.clabel(c3, inline=True, fontsize=8, fmt="%1.0f%%")
    ax3.set_xscale("log")
    ax3.set_yscale("log")
    ax3.set_title("c   F = 1.25 (Expansión tisular 25%): Penetrancia a 80 años", fontsize=11, fontweight="bold")
    ax3.set_xlabel("d (log scale)")
    ax3.set_ylabel(r"$\xi = \beta/\mu$ (log scale)")
    fig.colorbar(im3, ax=ax3)

    # F=1.25: Shift 50 vs 80
    ax4 = axes[1, 1]
    grid_125_50 = df_125.pivot(index="xi", columns="d", values="risk_50").values
    ax4.pcolormesh(d_vals, xi_vals, grid_125_80, shading="gouraud", cmap="Blues", vmin=0, vmax=100, alpha=0.8)
    cs80_2 = ax4.contour(d_vals, xi_vals, grid_125_80, levels=[50, 70, 90], colors="black", linestyles="solid")
    cs50_2 = ax4.contour(d_vals, xi_vals, grid_125_50, levels=[20, 46, 60], colors="crimson", linestyles="dashed")
    ax4.clabel(cs80_2, inline=True, fontsize=8, fmt="80y: %1.0f%%")
    ax4.clabel(cs50_2, inline=True, fontsize=8, fmt="50y: %1.0f%%")
    ax4.set_xscale("log")
    ax4.set_yscale("log")
    ax4.set_title("d   F = 1.25: Desplazamiento 80a (sólido) vs 50a (discontinuo)", fontsize=11, fontweight="bold")
    ax4.set_xlabel("d (log scale)")
    ax4.set_ylabel(r"$\xi = \beta/\mu$ (log scale)")

    plt.tight_layout()
    out_png = FIG_DIR / "Figure_grid_xi_d_compare_F.png"
    out_pdf = FIG_DIR / "Figure_grid_xi_d_compare_F.pdf"
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.savefig(out_pdf, bbox_inches="tight")
    plt.close()
    print(f"\n[DONE] Comparativa F=0.75 y F=1.25 guardada en: {out_png}")


if __name__ == "__main__":
    main()
