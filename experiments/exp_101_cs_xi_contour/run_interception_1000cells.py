#!/usr/bin/env python3
"""
run_interception_1000cells.py — Cancer interception experiment under 1000 cells neutral dynamics.
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
DOSES = (0.10, 0.25, 0.50)


def make_design(post: pd.DataFrame, out_path: Path, n_base_per_age: int = 1000, seed: int = 20260922) -> pd.DataFrame:
    print(f"Generating matched intervention design across ages {AGES} and doses {DOSES} (base per age: {n_base_per_age})...")
    rng = np.random.default_rng(seed)
    w = post["weight"].to_numpy(float, copy=True)
    w = w / w.sum()
    rows = []
    tissue_global = 0

    n_high_dose = int(n_base_per_age * 0.75)

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
                if dose > 0.10 and i >= n_high_dose:
                    continue
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
    print("Computing matched summaries, ARR, dose-response, and interaction contrasts...")

    # 1. Figure 5a: Risk curves for age 30, 50% attenuation
    curves = []
    for condition in ("none", "D1", "D2", "D1+D2"):
        if condition == "none":
            z = res[(res["intervention_age"] == 30) & (res["condition"] == "none") & (res["eligible"] == 1)]
            tids = res[(res["intervention_age"] == 30) & (res["condition"] == "D1") & np.isclose(res["dose"], 0.50)]["tissue_id"]
            z = z[z["tissue_id"].isin(tids)]
        else:
            z = res[(res["intervention_age"] == 30) & (res["condition"] == condition) & np.isclose(res["dose"], 0.50) & (res["eligible"] == 1)]
        c = cumulative_curve(z, 30, 80)
        c["condition"] = condition
        curves.append(c)
    curves_df = pd.concat(curves, ignore_index=True)
    curves_df.to_csv(out_dir / "Figure5a_risk_curves.csv", index=False)
    print("  -> Saved Figure5a_risk_curves.csv")

    # 2. Figure 5b: ARR by age 70 for 50% dose across start ages 30, 40, 50
    arr = []
    tests = []
    for age in AGES:
        for condition in ("D1", "D2", "D1+D2"):
            m = matched_set(res, age, 0.50, condition)
            r0 = 100.0 * m["tumor_by70_none"].mean()
            ri = 100.0 * m["tumor_by70_tr"].mean()
            arr.append({
                "intervention_age": age,
                "condition": condition,
                "dose": 0.50,
                "risk_none_70": r0,
                "risk_intervention_70": ri,
                "ARR_pp": r0 - ri,
                "n": len(m)
            })
            table = pd.crosstab(m["tumor_by70_none"], m["tumor_by70_tr"]).reindex(index=[0, 1], columns=[0, 1], fill_value=0)
            test = mcnemar(table.to_numpy(), exact=False, correction=False)
            tests.append({
                "analysis": "McNemar_age70",
                "intervention_age": age,
                "condition": condition,
                "dose": 0.50,
                "statistic": float(test.statistic),
                "p": float(test.pvalue),
                "n": len(m)
            })
    arr_df = pd.DataFrame(arr)
    arr_df.to_csv(out_dir / "Figure5b_ARR_age70_50pct.csv", index=False)
    pd.DataFrame(tests).to_csv(out_dir / "Figure5b_McNemar_tests.csv", index=False)
    print("  -> Saved Figure5b_ARR_age70_50pct.csv and Figure5b_McNemar_tests.csv")

    # 3. Figure 5c: Dose-response at age 30
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
    pd.DataFrame(dose_rows).to_csv(out_dir / "Figure5c_dose_response_age30.csv", index=False)
    pd.DataFrame(gee_rows).to_csv(out_dir / "Figure5c_GEE_tests.csv", index=False)
    print("  -> Saved Figure5c_dose_response_age30.csv and Figure5c_GEE_tests.csv")

    # 4. Figure S5: Interaction on absolute-risk scale through age 80
    rng = np.random.default_rng(seed)
    inter = []
    for age in AGES:
        for dose in DOSES:
            d1 = matched_set(res, age, dose, "D1")
            d2 = matched_set(res, age, dose, "D2")
            both = matched_set(res, age, dose, "D1+D2")
            x = d1[["tissue_id", "tumor_by80_none", "tumor_by80_tr"]].rename(columns={"tumor_by80_tr": "y_d1"})
            x = x.merge(d2[["tissue_id", "tumor_by80_tr"]].rename(columns={"tumor_by80_tr": "y_d2"}), on="tissue_id")
            x = x.merge(both[["tissue_id", "tumor_by80_tr"]].rename(columns={"tumor_by80_tr": "y_both"}), on="tissue_id")
            y0 = x["tumor_by80_none"].to_numpy(float)
            y1 = x["y_d1"].to_numpy(float)
            y2 = x["y_d2"].to_numpy(float)
            yb = x["y_both"].to_numpy(float)

            def icontrast(ii):
                r0 = y0[ii].mean()
                a1 = r0 - y1[ii].mean()
                a2 = r0 - y2[ii].mean()
                ab = r0 - yb[ii].mean()
                return 100.0 * (ab - a1 - a2)

            base = np.arange(len(x))
            est = icontrast(base)
            bs = np.empty(boot)
            for b in range(boot):
                ii = rng.integers(0, len(x), len(x))
                bs[b] = icontrast(ii)
            inter.append({
                "intervention_age": age,
                "dose": dose,
                "interaction_pp": est,
                "ci_lo": float(np.quantile(bs, 0.025)),
                "ci_hi": float(np.quantile(bs, 0.975)),
                "n": len(x)
            })
    pd.DataFrame(inter).to_csv(out_dir / "FigureS5_interaction.csv", index=False)
    print("  -> Saved FigureS5_interaction.csv")


def make_design_75_only(existing_design_path: Path, out_path_75: Path) -> pd.DataFrame:
    print(f"Reading base tissues from {existing_design_path}...")
    df_prev = pd.read_csv(existing_design_path)
    base_tissues = df_prev[df_prev["condition"] == "none"].copy()

    rows = []
    for _, b in base_tissues.iterrows():
        base_dict = b.to_dict()
        base_dict["dose"] = 0.75
        base_dict["r_d1"] = 0.75
        base_dict["r_d2"] = 0.0
        base_dict["condition"] = "D1"
        rows.append(dict(base_dict))

        base_dict["r_d1"] = 0.0
        base_dict["r_d2"] = 0.75
        base_dict["condition"] = "D2"
        rows.append(dict(base_dict))

        base_dict["r_d1"] = 0.75
        base_dict["r_d2"] = 0.75
        base_dict["condition"] = "D1+D2"
        rows.append(dict(base_dict))

    df_75 = pd.DataFrame(rows)
    df_75.to_csv(out_path_75, index=False)
    print(f"Design 75% generated: exactly {len(df_75):,} conditions saved to {out_path_75}")
    return df_75


def summarize_with_75(res_merged: pd.DataFrame, out_dir: Path):
    print("Computing combined summaries for Figure 5 (a, b, c, d)...")
    ALL_DOSES = (0.10, 0.25, 0.50, 0.75)

    # 1. Figure 5a: Risk curves for age 30, 50% dose
    curves_50 = []
    for condition in ("none", "D1", "D2", "D1+D2"):
        if condition == "none":
            z = res_merged[(res_merged["intervention_age"] == 30) & (res_merged["condition"] == "none") & (res_merged["eligible"] == 1)]
            tids = res_merged[(res_merged["intervention_age"] == 30) & (res_merged["condition"] == "D1") & np.isclose(res_merged["dose"], 0.50)]["tissue_id"]
            z = z[z["tissue_id"].isin(tids)]
        else:
            z = res_merged[(res_merged["intervention_age"] == 30) & (res_merged["condition"] == condition) & np.isclose(res_merged["dose"], 0.50) & (res_merged["eligible"] == 1)]
        c = cumulative_curve(z, 30, 80)
        c["condition"] = condition
        curves_50.append(c)
    pd.concat(curves_50, ignore_index=True).to_csv(out_dir / "Figure5a_risk_curves_50pct.csv", index=False)

    # 2. Figure 5b: Risk curves for age 30, 75% dose
    curves_75 = []
    for condition in ("none", "D1", "D2", "D1+D2"):
        if condition == "none":
            z = res_merged[(res_merged["intervention_age"] == 30) & (res_merged["condition"] == "none") & (res_merged["eligible"] == 1)]
            tids = res_merged[(res_merged["intervention_age"] == 30) & (res_merged["condition"] == "D1") & np.isclose(res_merged["dose"], 0.75)]["tissue_id"]
            z = z[z["tissue_id"].isin(tids)]
        else:
            z = res_merged[(res_merged["intervention_age"] == 30) & (res_merged["condition"] == condition) & np.isclose(res_merged["dose"], 0.75) & (res_merged["eligible"] == 1)]
        c = cumulative_curve(z, 30, 80)
        c["condition"] = condition
        curves_75.append(c)
    pd.concat(curves_75, ignore_index=True).to_csv(out_dir / "Figure5b_risk_curves_75pct.csv", index=False)

    # 3. Figure 5c: ARR by age 70 for 75% dose across start ages 30, 40, 50
    arr_75 = []
    tests_75 = []
    for age in AGES:
        for condition in ("D1", "D2", "D1+D2"):
            m = matched_set(res_merged, age, 0.75, condition)
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
                "analysis": "McNemar_age70_75pct",
                "intervention_age": age,
                "condition": condition,
                "dose": 0.75,
                "statistic": float(test.statistic),
                "p": float(test.pvalue),
                "n": len(m)
            })
    pd.DataFrame(arr_75).to_csv(out_dir / "Figure5c_ARR_age70_75pct.csv", index=False)
    pd.DataFrame(tests_75).to_csv(out_dir / "Figure5c_McNemar_tests_75pct.csv", index=False)

    # 4. Figure 5d: Dose-response at age 30 across (10%, 25%, 50%, 75%)
    dose_rows = []
    gee_rows = []
    for condition in ("D1", "D2", "D1+D2"):
        pieces = []
        for dose in ALL_DOSES:
            m = matched_set(res_merged, 30, dose, condition)
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
            print(f"GEE warning: {e}")
    pd.DataFrame(dose_rows).to_csv(out_dir / "Figure5d_dose_response_age30.csv", index=False)
    pd.DataFrame(gee_rows).to_csv(out_dir / "Figure5d_GEE_tests.csv", index=False)
    print("  -> Combined summaries successfully written.")


def main() -> int:
    parser = argparse.ArgumentParser(description="Cancer interception runner for 1000 cells.")
    parser.add_argument("--workers", type=int, default=12, help="OMP worker threads (default: 12)")
    parser.add_argument("--base-tissues", type=int, default=1000, help="Base tissues per start age (default: 1000)")
    parser.add_argument("--cells", type=int, default=1000, help="Initial cells per tissue (default: 1000)")
    parser.add_argument("--binary", type=str, default="./build/run_paper_interventions")
    parser.add_argument("--post-csv", type=str, default="experiments/exp_101_cs_xi_contour/runs/run_abc_xi_wide_20261120_171832/particles_gen_10.csv")
    parser.add_argument("--output-dir", type=str, default="experiments/exp_101_cs_xi_contour/results_1000cells_final/interception")
    parser.add_argument("--only-75", action="store_true", help="Simulate only 75% conditions and merge with existing results")
    args = parser.parse_args()

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    if args.only_75:
        prev_design = out_dir / "interception_design.csv"
        prev_results = out_dir / "interception_results.csv"
        design_75 = out_dir / "interception_design_75.csv"
        results_75 = out_dir / "interception_results_75.csv"
        merged_results = out_dir / "interception_results_merged.csv"

        make_design_75_only(prev_design, design_75)
        run_binary(args.binary, design_75, results_75, cells=args.cells, threads=args.workers)

        print("Combinando outcomes...")
        df_prev_res = pd.read_csv(prev_results)
        df_75_res = pd.read_csv(results_75)
        df_merged = pd.concat([df_prev_res, df_75_res], ignore_index=True)
        df_merged.to_csv(merged_results, index=False)
        print(f"Total combinado: {len(df_merged):,} outcomes.")
        summarize_with_75(df_merged, out_dir)
        return 0

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
