#!/usr/bin/env python3
"""
simulate_victor_8cases.py — Evaluates the 8 parameter configurations proposed by Víctor.
Runs 200 Monte Carlo tissues for penetrance (Kuchenbaecker comparison) and
100 virtual tissues for stage dynamics (testing 75-100% baseline retention).

Configurations:
  Case A: mu=0.003,  d=0.05,   beta in [0.04475, 0.04818] (mid=0.046465)
  Case B: mu=0.003,  d=0.02,   beta in [0.01666, 0.02019] (mid=0.018425)
  Case C: mu=0.003,  d=0.01,   beta in [0.00692, 0.01049] (mid=0.008705)
  Case D: mu=0.003,  d=0.006,  beta in [0.00297, 0.00655] (mid=0.004760)
  Case E: mu=0.0015, d=0.006,  beta in [0.00447, 0.00804] (mid=0.006255)
  Case F: mu=0.0015, d=0.0045, beta in [0.00298, 0.00656] (mid=0.004770)
  Case G: mu=0.001,  d=0.004,  beta in [0.00299, 0.00657] (mid=0.004780)
  Case H: mu=0.001,  d=0.003,  beta in [0.00199, 0.00558] (mid=0.003785)
"""

from __future__ import annotations

import os
import re
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

plt.switch_backend("Agg")

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parents[2]
DATA_DIR = SCRIPT_DIR / "data"
FIG_DIR = SCRIPT_DIR / "figures"

BOOTSTRAP_BIN = PROJECT_ROOT / "build" / "run_bootstrapping"
STAGE_BIN = PROJECT_ROOT / "build" / "run_stage_trajectories"

CASES = [
    {"label": "A", "mu": 0.0030, "d": 0.0500, "b_lo": 0.04475, "b_hi": 0.04818, "b_mid": 0.046465},
    {"label": "B", "mu": 0.0030, "d": 0.0200, "b_lo": 0.01666, "b_hi": 0.02019, "b_mid": 0.018425},
    {"label": "C", "mu": 0.0030, "d": 0.0100, "b_lo": 0.00692, "b_hi": 0.01049, "b_mid": 0.008705},
    {"label": "D", "mu": 0.0030, "d": 0.0060, "b_lo": 0.00297, "b_hi": 0.00655, "b_mid": 0.004760},
    {"label": "E", "mu": 0.0015, "d": 0.0060, "b_lo": 0.00447, "b_hi": 0.00804, "b_mid": 0.006255},
    {"label": "F", "mu": 0.0015, "d": 0.0045, "b_lo": 0.00298, "b_hi": 0.00656, "b_mid": 0.004770},
    {"label": "G", "mu": 0.0010, "d": 0.0040, "b_lo": 0.00299, "b_hi": 0.00657, "b_mid": 0.004780},
    {"label": "H", "mu": 0.0010, "d": 0.0030, "b_lo": 0.00199, "b_hi": 0.00558, "b_mid": 0.003785},
]

KUCHENBAECKER = {
    30: (4.0, (2.0, 8.0)),
    40: (24.0, (20.0, 29.0)),
    50: (43.0, (38.0, 48.0)),
    60: (56.0, (51.0, 61.0)),
    70: (65.0, (60.0, 70.0)),
    80: (72.0, (65.0, 79.0)),
}


def run_penetrance(case: dict, n_runs: int = 300) -> dict:
    """Runs run_bootstrapping and parses risk at milestone ages."""
    cmd = [
        str(BOOTSTRAP_BIN),
        "--n-runs", str(n_runs),
        "--weighted-sse",
        "--big-bang",
        "--n-cells", "500",
        "--max-healthy-cells", "500",
        "--division-rate", str(case["d"]),
        "--tp53-rate", str(case["mu"]),
        "--brca1-rate", str(case["b_mid"]),
        "--low-delta", "0.102",
        "--high-delta", "0.200",
        "--neoplastic-div-rate", "0.250",
        "--output", "/dev/null",
        "--seed-offset", "42000",
    ]
    env = dict(os.environ)
    env["OMP_NUM_THREADS"] = "1"
    res = subprocess.run(cmd, capture_output=True, text=True, env=env)

    dist_re = re.search(r"weighted SSE\s*=\s*([\d.eE+\-]+)", res.stdout)
    sse_val = float(dist_re.group(1)) if dist_re else float("nan")

    risks = {}
    for m in re.finditer(r"^\s*(\d{2})\s+([\d.]+)\s+", res.stdout, re.MULTILINE):
        risks[int(m.group(1))] = float(m.group(2))

    return {"case": case["label"], "sse": sse_val, "risks": risks}


def run_stage_trajectory(case: dict, n_tissues: int = 100) -> Path:
    """Runs run_stage_trajectories and saves CSV."""
    csv_path = DATA_DIR / f"stage_trajectories_case_{case['label']}.csv"
    cmd = [
        str(STAGE_BIN),
        "--n_tissues", str(n_tissues),
        "--n_cells", "500",
        "--max_cells", "500",
        "--max_t", "80",
        "--division_rate", str(case["d"]),
        "--brca1_rate", str(case["b_mid"]),
        "--tp53_rate", str(case["mu"]),
        "--delta_low", "0.102",
        "--delta_high", "0.200",
        "--neoplastic_div_rate", "0.250",
        "--output", str(csv_path),
        "--seed", "42",
    ]
    subprocess.run(cmd, capture_output=True, text=True, check=True)
    return csv_path


def main():
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    print("=" * 70)
    print("Evaluating Víctor's 8 Parameter Configurations (A to H)")
    print("=" * 70)

    # 1. Run Penetrance Simulations in parallel
    print("\n[Step 1/3] Running penetrance evaluations (N=300 runs each)...")
    with ThreadPoolExecutor(max_workers=8) as executor:
        penetrance_results = list(executor.map(run_penetrance, CASES))

    df_pen = []
    for r in penetrance_results:
        row = {"Case": r["case"], "Weighted_SSE": r["sse"]}
        for age in sorted(KUCHENBAECKER.keys()):
            row[f"Risk_{age}"] = r["risks"].get(age, np.nan)
        df_pen.append(row)
    df_pen = pd.DataFrame(df_pen)
    df_pen.to_csv(DATA_DIR / "victor_8cases_penetrance_summary.csv", index=False)
    print("\nPenetrance Fit Summary:")
    print(df_pen.to_string(index=False))

    # 2. Run Stage Dynamics Simulations in parallel
    print("\n[Step 2/3] Simulating Stage Trajectories (100 tissues each)...")
    with ThreadPoolExecutor(max_workers=8) as executor:
        stage_files = list(executor.map(run_stage_trajectory, CASES))

    # 3. Render Combined Comparison Figure
    print("\n[Step 3/3] Generating comparative multi-panel figure...")
    fig, axes = plt.subplots(2, 2, figsize=(14, 11), dpi=300)

    colors = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f"]
    markers = ["o", "s", "^", "v", "D", "P", "X", "*"]

    # Panel (a): Penetrance curves vs Kuchenbaecker
    ax_pen = axes[0, 0]
    ages = np.array(sorted(KUCHENBAECKER.keys()))
    kuch_vals = [KUCHENBAECKER[a][0] for a in ages]
    kuch_lo = [KUCHENBAECKER[a][1][0] for a in ages]
    kuch_hi = [KUCHENBAECKER[a][1][1] for a in ages]

    ax_pen.fill_between(ages, kuch_lo, kuch_hi, color="#9CA3AF", alpha=0.3, label="Kuchenbaecker 95% CI")
    ax_pen.plot(ages, kuch_vals, "k--", lw=2, marker="s", markersize=5, label="Observed (JAMA 2017)")

    for c_idx, c in enumerate(CASES):
        lbl = c["label"]
        row = df_pen[df_pen["Case"] == lbl].iloc[0]
        y_vals = [row[f"Risk_{a}"] for a in ages]
        sse = row["Weighted_SSE"]
        ax_pen.plot(ages, y_vals, color=colors[c_idx], marker=markers[c_idx], lw=1.6, markersize=5,
                    label=f"Case {lbl} (d={c['d']}, μ={c['mu']}, SSE={sse:.1f})")

    ax_pen.set_title("a   Cumulative Penetrance vs Kuchenbaecker 2017", fontsize=11, fontweight="bold", loc="left")
    ax_pen.set_xlabel("Age (years)", fontsize=10)
    ax_pen.set_ylabel("Cumulative Incidence (%)", fontsize=10)
    ax_pen.set_xlim(25, 85)
    ax_pen.set_ylim(0, 105)
    ax_pen.grid(True, linestyle="--", alpha=0.5)
    ax_pen.legend(fontsize=7.5, frameon=True, loc="upper left")

    # Panel (b): Normalized Baseline Healthy Cell Retention E[B(t)]/B(0)
    ax_base = axes[0, 1]
    ax_base.axhspan(0.75, 1.0, color="#DCFCE7", alpha=0.6, label="Víctor's Target Window (75–100%)")
    ax_base.axhline(0.75, color="#16A34A", linestyle=":", lw=1.2)
    ax_base.axhline(1.0, color="#16A34A", linestyle=":", lw=1.2)

    for c_idx, c in enumerate(CASES):
        lbl = c["label"]
        csv_p = DATA_DIR / f"stage_trajectories_case_{lbl}.csv"
        df_st = pd.read_csv(csv_p)
        b_mean = df_st.groupby("year")["baseline"].mean()
        b_ret = b_mean / 500.0
        ax_base.plot(b_mean.index, b_ret, color=colors[c_idx], lw=1.8, label=f"Case {lbl} (80y: {b_ret.loc[80]*100:.1f}%)")

    ax_base.set_title("b   Normalized Baseline Healthy Cell Retention $E[B(t)]/B(0)$", fontsize=11, fontweight="bold", loc="left")
    ax_base.set_xlabel("Tissue Age (years)", fontsize=10)
    ax_base.set_ylabel("Fraction of Initial Baseline Cells", fontsize=10)
    ax_base.set_xlim(1, 80)
    ax_base.set_ylim(0.0, 1.05)
    ax_base.grid(True, linestyle="--", alpha=0.5)
    ax_base.legend(fontsize=7.5, frameon=True, loc="lower left")

    # Panel (c): Mean Absolute Cell Counts (Baseline vs Tumoral) for 4 Selected Cases (A, B, E, H)
    ax_dyn = axes[1, 0]
    select_cases = ["A", "B", "E", "H"]
    for sc in select_cases:
        csv_p = DATA_DIR / f"stage_trajectories_case_{sc}.csv"
        df_st = pd.read_csv(csv_p)
        b_mean = df_st.groupby("year")["baseline"].mean()
        t_mean = df_st.groupby("year")["tumoral"].mean()
        c_case = next(c for c in CASES if c["label"] == sc)
        ax_dyn.plot(b_mean.index, b_mean, lw=1.8, label=f"Case {sc} Baseline (d={c_case['d']})")
        ax_dyn.plot(t_mean.index, t_mean, lw=1.8, linestyle="--", label=f"Case {sc} Tumoral")

    ax_dyn.set_title("c   Absolute Cell Counts: Baseline (solid) vs Tumoral (dashed)", fontsize=11, fontweight="bold", loc="left")
    ax_dyn.set_xlabel("Tissue Age (years)", fontsize=10)
    ax_dyn.set_ylabel("Mean Cells per Tissue (N0=500)", fontsize=10)
    ax_dyn.set_xlim(1, 80)
    ax_dyn.set_ylim(0, 550)
    ax_dyn.grid(True, linestyle="--", alpha=0.5)
    ax_dyn.legend(fontsize=7.5, frameon=True, loc="upper right")

    # Panel (d): Comparative Parameter Summary Table
    ax_tbl = axes[1, 1]
    ax_tbl.axis("off")
    ax_tbl.set_title("d   Víctor's Parameter Table & Penetrance Error (Weighted SSE)", fontsize=11, fontweight="bold", loc="left")

    table_data = []
    for c in CASES:
        lbl = c["label"]
        row = df_pen[df_pen["Case"] == lbl].iloc[0]
        csv_p = DATA_DIR / f"stage_trajectories_case_{lbl}.csv"
        df_st = pd.read_csv(csv_p)
        ret_80 = (df_st.groupby("year")["baseline"].mean().loc[80] / 500.0) * 100.0
        table_data.append([
            lbl,
            f"{c['mu']:.4f}",
            f"{c['d']:.4f}",
            f"{c['b_mid']:.5f}",
            f"[{c['b_lo']:.4f}, {c['b_hi']:.4f}]",
            f"{ret_80:.1f}%",
            f"{row['Weighted_SSE']:.1f}",
        ])

    columns = ["Case", "μ (TP53)", "d (div)", "β (mid)", "β range", "80y Ret.", "Weighted SSE"]
    table = ax_tbl.table(cellText=table_data, colLabels=columns, loc="center", cellLoc="center")
    table.auto_set_font_size(False)
    table.set_fontsize(8.5)
    table.scale(1.05, 1.6)

    # Style table header
    for (row_idx, col_idx), cell in table.get_celld().items():
        if row_idx == 0:
            cell.set_facecolor("#E5E7EB")
            cell.set_text_props(fontweight="bold")
        elif col_idx == 0:
            cell.set_text_props(fontweight="bold")

    plt.tight_layout()
    out_png = FIG_DIR / "Figure_victor_8cases_comparison.png"
    out_pdf = FIG_DIR / "Figure_victor_8cases_comparison.pdf"
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.savefig(out_pdf, bbox_inches="tight")
    plt.close()

    print(f"\n[DONE] Comparative analysis completed successfully!")
    print(f"  -> PNG: {out_png}")
    print(f"  -> PDF: {out_pdf}")


if __name__ == "__main__":
    main()
