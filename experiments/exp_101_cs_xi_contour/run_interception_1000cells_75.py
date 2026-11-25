#!/usr/bin/env python3
"""
run_interception_1000cells_75.py — Cancer interception experiment including 75% attenuation.
Doses: 10%, 25%, 50%, 75% across ages 30, 40, 50.
"""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from statsmodels.genmod.cov_struct import Exchangeable
from statsmodels.stats.contingency_tables import mcnemar

AGES = (30, 40, 50)
DOSES = (0.10, 0.25, 0.50, 0.75)


def make_design(post: pd.DataFrame, out_path: Path, n_base_per_age: int = 1000, seed: int = 20260922) -> pd.DataFrame:
    print(f"Generating matched intervention design across ages {AGES} and doses {DOSES} (base per age: {n_base_per_age})...")
    rng = np.random.default_rng(seed)
    w = post["weight"].to_numpy(float, copy=True)
    w = w / w.sum()
    rows = []
    tissue_global = 0

    for age in AGES:
        idx = rng.choice(len(post), size=n_base_per_age, replace=True, p=w)
        for i, pidx in enumerate(idx):
            r = post.iloc[pidx]
            tid = tissue_global + i
            tissue_seed = 10000 + (age // 10) * 10000 + i
            base = {
                "tissue_id": tid,
                "seed": tissue_seed,
                "brca1_rate": float(r["beta"]),
                "tp53_rate": float(r["mu"]),
                "division_rate": float(r["division_rate"]),
                "low_delta": float(r["low_delta"]),
                "high_delta": float(r["high_delta"]),
                "neoplastic_div_rate": float(r["neoplastic_div_rate"]),
                "d1_threshold": 2.0,
                "d2_threshold": 5.0,
                "tumor_threshold": 0.05,
                "kappa_brca1": 2.0,
                "age_drift": 1e-5,
                "intervention_age": age,
            }
            rows.append({**base, "condition": "none", "dose": 0.0, "r_d1": 0.0, "r_d2": 0.0})
            for dose in DOSES:
                rows.append({**base, "condition": "D1", "dose": dose, "r_d1": dose, "r_d2": 0.0})
                rows.append({**base, "condition": "D2", "dose": dose, "r_d1": 0.0, "r_d2": dose})
                rows.append({**base, "condition": "D1+D2", "dose": dose, "r_d1": dose, "r_d2": dose})
        tissue_global += 10000

    df = pd.DataFrame(rows)
    df.to_csv(out_path, index=False)
    print(f"Design matrix generated: {len(df):,} combinations saved to {out_path}")
    return df


def run_binary(binary_path: str, design_path: Path, result_path: Path, cells: int = 1000, threads: int = 12):
    print(f"Launching C++ binary {binary_path} with OMP_NUM_THREADS={threads}, cells={cells}...")
    env = dict(os.environ)
    env["OMP_NUM_THREADS"] = str(threads)
    t0 = time.time()
    cmd = [binary_path, "--design", str(design_path), "--output", str(result_path), "--cells", str(cells)]
    subprocess.check_call(cmd, env=env)
    elapsed = time.time() - t0
    print(f"[SUCCESS] C++ simulation finished in {elapsed:.2f}s! Outputs: {result_path}")


def cumulative_curve(d: pd.DataFrame, age0: int, end: int = 80) -> pd.DataFrame:
    ages = np.arange(age0, end + 1)
    return pd.DataFrame({
        "age": ages,
        "risk": [100.0 * np.mean((d["onset_year"] >= 0) & (d["onset_year"] <= a)) for a in ages]
    })


def matched_set(res: pd.DataFrame, age: int, dose: float, condition: str) -> pd.DataFrame:
    none = res[(res["intervention_age"] == age) & (res["condition"] == "none")][
        ["tissue_id", "eligible", "tumor_by70", "tumor_by80", "onset_year"]
    ]
    tr = res[(res["intervention_age"] == age) & (res["condition"] == condition) & np.isclose(res["dose"], dose)][
        ["tissue_id", "eligible", "tumor_by70", "tumor_by80", "onset_year"]
    ]
    m = none.merge(tr, on="tissue_id", suffixes=("_none", "_tr"))
    if not np.array_equal(m["eligible_none"].to_numpy(), m["eligible_tr"].to_numpy()):
        raise RuntimeError("Eligibility mismatch between matched conditions.")
    return m[m["eligible_none"] == 1].copy()


def summarize(res: pd.DataFrame, out_dir: Path, boot: int = 500, seed: int = 123):
    print("Computing matched summaries for 50% and 75% doses...")

    # 1. Figure 5a: Risk curves for age 30, 50% attenuation
    curves_50 = []
    for condition in ("none", "D1", "D2", "D1+D2"):
        if condition == "none":
            z = res[(res["intervention_age"] == 30) & (res["condition"] == "none") & (res["eligible"] == 1)]
            tids = res[(res["intervention_age"] == 30) & (res["condition"] == "D1") & np.isclose(res["dose"], 0.50)]["tissue_id"]
            z = z[z["tissue_id"].isin(tids)]
        else:
            z = res[(res["intervention_age"] == 30) & (res["condition"] == condition) & np.isclose(res["dose"], 0.50) & (res["eligible"] == 1)]
        c = cumulative_curve(z, 30, 80)
        c["condition"] = condition
        curves_50.append(c)
    curves_50_df = pd.concat(curves_50, ignore_index=True)
    curves_50_df.to_csv(out_dir / "Figure5a_risk_curves_50pct.csv", index=False)
    print("  -> Saved Figure5a_risk_curves_50pct.csv")

    # 2. Figure 5b: Risk curves for age 30, 75% attenuation
    curves_75 = []
    for condition in ("none", "D1", "D2", "D1+D2"):
        if condition == "none":
            z = res[(res["intervention_age"] == 30) & (res["condition"] == "none") & (res["eligible"] == 1)]
            tids = res[(res["intervention_age"] == 30) & (res["condition"] == "D1") & np.isclose(res["dose"], 0.75)]["tissue_id"]
            z = z[z["tissue_id"].isin(tids)]
        else:
            z = res[(res["intervention_age"] == 30) & (res["condition"] == condition) & np.isclose(res["dose"], 0.75) & (res["eligible"] == 1)]
        c = cumulative_curve(z, 30, 80)
        c["condition"] = condition
        curves_75.append(c)
    curves_75_df = pd.concat(curves_75, ignore_index=True)
    curves_75_df.to_csv(out_dir / "Figure5b_risk_curves_75pct.csv", index=False)
    print("  -> Saved Figure5b_risk_curves_75pct.csv")

    # 3. Figure 5c: Effect of intervention timing (ARR by age 70 for 75% dose across ages 30, 40, 50)
    arr_75 = []
    tests_75 = []
    for age in AGES:
        for condition in ("D1", "D2", "D1+D2"):
            m = matched_set(res, age, 0.75, condition)
            r0 = 100.0 * m["tumor_by70_none"].mean()
            ri = 100.0 * m["tumor_by70_tr"].mean()
            arr_75.append({
                "intervention_age": age,
                "condition": condition,
                "dose": 0.75,
                "risk_none_70": r0,
                "risk_intervention_70": ri,
                "ARR_pp": r0 - ri,
                "n": len(m)
            })
            table = pd.crosstab(m["tumor_by70_none"], m["tumor_by70_tr"]).reindex(index=[0, 1], columns=[0, 1], fill_value=0)
            test = mcnemar(table.to_numpy(), exact=False, correction=False)
            tests_75.append({
                "analysis": "McNemar_age70",
                "intervention_age": age,
                "condition": condition,
                "dose": 0.75,
                "statistic": float(test.statistic),
                "p": float(test.pvalue),
                "n": len(m)
            })
    pd.DataFrame(arr_75).to_csv(out_dir / "Figure5c_ARR_age70_75pct.csv", index=False)
    pd.DataFrame(tests_75).to_csv(out_dir / "Figure5c_McNemar_tests_75pct.csv", index=False)
    print("  -> Saved Figure5c_ARR_age70_75pct.csv and Figure5c_McNemar_tests_75pct.csv")

    # 4. Figure 5d: Dose-response at age 30 across (10%, 25%, 50%, 75%)
    dose_rows = []
    gee_rows = []
    for condition in ("D1", "D2", "D1+D2"):
        pieces = []
        for dose in DOSES:
            m = matched_set(res, 30, dose, condition)
            r0 = 100.0 * m["tumor_by70_none"].mean()
            ri = 100.0 * m["tumor_by70_tr"].mean()
            dose_rows.append({
                "condition": condition,
                "dose": dose,
                "risk_none_70": r0,
                "risk_intervention_70": ri,
                "ARR_pp": r0 - ri,
                "n": len(m)
            })
            x = m[["tissue_id", "tumor_by70_tr"]].copy()
            x["dose"] = dose
            x.rename(columns={"tumor_by70_tr": "tumor"}, inplace=True)
            pieces.append(x)
        g = pd.concat(pieces, ignore_index=True)
        try:
            model = smf.gee("tumor ~ dose", groups="tissue_id", data=g,
                            family=sm.families.Binomial(), cov_struct=Exchangeable()).fit()
            gee_rows.append({
                "condition": condition,
                "coef_dose": float(model.params["dose"]),
                "se": float(model.bse["dose"]),
                "p": float(model.pvalues["dose"])
            })
        except Exception as e:
            print(f"Warning GEE: {e}")
    pd.DataFrame(dose_rows).to_csv(out_dir / "Figure5d_dose_response_age30.csv", index=False)
    pd.DataFrame(gee_rows).to_csv(out_dir / "Figure5d_GEE_tests.csv", index=False)
    print("  -> Saved Figure5d_dose_response_age30.csv and Figure5d_GEE_tests.csv")


def main() -> int:
    parser = argparse.ArgumentParser(description="Cancer interception runner including 75% dose.")
    parser.add_argument("--workers", type=int, default=12, help="OMP worker threads (default: 12)")
    parser.add_argument("--base-tissues", type=int, default=1000, help="Base tissues per start age (default: 1000)")
    parser.add_argument("--cells", type=int, default=1000, help="Initial cells per tissue (default: 1000)")
    parser.add_argument("--binary", type=str, default="./build/run_paper_interventions")
    parser.add_argument("--post-csv", type=str, default="experiments/exp_101_cs_xi_contour/runs/run_abc_xi_wide_20261120_171832/particles_gen_10.csv")
    parser.add_argument("--output-dir", type=str, default="experiments/exp_101_cs_xi_contour/results_1000cells_final/interception_75")
    args = parser.parse_args()

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    post = pd.read_csv(args.post_csv)
    design_path = out_dir / "interception_design.csv"
    results_path = out_dir / "interception_results.csv"

    make_design(post, design_path, n_base_per_age=args.base_tissues)
    run_binary(args.binary, design_path, results_path, cells=args.cells, threads=args.workers)
    res = pd.read_csv(results_path)
    summarize(res, out_dir)
    print(f"\n[DONE] Interception simulation and summary completed in: {out_dir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
