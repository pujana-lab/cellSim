#!/usr/bin/env python3
"""
run_homeostasis_features.py — Evaluates emergent tumor features across posterior particles
from the homeostatic ABC-SMC calibration (1000 cells, p in [0.047, 0.057], nu = 10*p).
Simulates virtual tissues using build/run_paper_features.
"""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
import tempfile
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Tuple

import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency, kendalltau

AGE_ORDER = ["<30", "30-39", "40-49", "50-59", "60-69", "70-79", "80"]
GENO_ORDER = ["Exclusively BRCA1+/-", "Mixed", "Exclusively BRCA1-/-"]
FOUNDER_ORDER = ["1", "2", "3", ">=4"]


def age_group(onset: int) -> str:
    if onset < 30:
        return "<30"
    if onset < 40:
        return "30-39"
    if onset < 50:
        return "40-49"
    if onset < 60:
        return "50-59"
    if onset < 70:
        return "60-69"
    if onset < 80:
        return "70-79"
    return "80"


def run_one_particle(item: Tuple[int, pd.Series], binary: str, n_per: int) -> pd.DataFrame:
    i, row = item
    with tempfile.NamedTemporaryFile(suffix=".csv", delete=False) as t:
        temp_path = t.name

    div_rate = float(row["division_rate"])
    neop_rate = float(row.get("neoplastic_div_rate", 10.0 * div_rate))

    cmd = [
        binary,
        "--n-runs", str(n_per),
        "--seed-offset", str(20000 + i * n_per),
        "--n-cells", "600",
        "--max-healthy-cells", "600",
        "--brca1-rate", str(float(row["brca1_rate"])),
        "--low-delta", str(float(row["low_delta"])),
        "--high-delta", str(float(row["high_delta"])),
        "--division-rate", str(div_rate),
        "--neoplastic-div-rate", str(neop_rate),
        "--tp53-rate", "0.003",
        "--d1-threshold", "2.0",
        "--d2-threshold", "5.0",
        "--tumor-threshold", "0.05",
        "--kappa-brca1", "2.0",
        "--age-drift", "0.00001",
        "--output", temp_path,
    ]
    env = dict(os.environ)
    env["OMP_NUM_THREADS"] = "1"
    res = subprocess.run(cmd, capture_output=True, text=True, env=env)
    if res.returncode != 0:
        if os.path.exists(temp_path):
            os.unlink(temp_path)
        raise RuntimeError(f"Error in run_paper_features: {res.stderr}")

    df = pd.read_csv(temp_path)
    if os.path.exists(temp_path):
        os.unlink(temp_path)

    df["particle"] = i
    w = float(row.get("weight", 1.0))
    df["posterior_weight"] = w
    df["tissue_weight"] = w / n_per
    return df


def simulate_cohort(post_csv: Path, binary: str, n_per: int, workers: int, out_dir: Path) -> pd.DataFrame:
    post = pd.read_csv(post_csv)
    if "weight" in post.columns and post["weight"].sum() > 0:
        post["weight"] = post["weight"] / post["weight"].sum()
    else:
        post["weight"] = 1.0 / len(post)

    print(f"=========================================================================")
    print(f"Simulating Emergent Features: {len(post)} particles x {n_per} tissues = {len(post)*n_per:,} tissues")
    print(f"Binary: {binary} | Workers: {workers}")
    print(f"=========================================================================")

    t0 = time.time()
    jobs = list(post.iterrows())
    dfs = []
    with ThreadPoolExecutor(max_workers=workers) as ex:
        for k, d in enumerate(ex.map(lambda item: run_one_particle(item, binary, n_per), jobs), 1):
            dfs.append(d)
            if k % 25 == 0 or k == len(jobs):
                print(f"  -> Processed {k}/{len(jobs)} particles ({k*n_per:,} tissues)...", flush=True)

    elapsed = time.time() - t0
    all_df = pd.concat(dfs, ignore_index=True)
    raw_csv = out_dir / f"emergent_{len(post)*n_per}_tissues.csv"
    all_df.to_csv(raw_csv, index=False)
    print(f"[SUCCESS] Simulation completed in {elapsed:.2f}s! Saved: {raw_csv}")
    return all_df


def summarize_features(raw_df: pd.DataFrame, out_dir: Path):
    d = raw_df[raw_df["onset_year"] >= 0].copy()
    d["age_group"] = d["onset_year"].map(age_group)
    d["founder_category"] = d["founder_count"].map(lambda x: ">=4" if x >= 4 else str(int(x)))

    mapping = {
        "BRCA1+/-": "Exclusively BRCA1+/-",
        "mixed": "Mixed",
        "BRCA1-/-": "Exclusively BRCA1-/-",
    }
    d["genotype_cat"] = d["brca1_category"].map(mapping)

    total_w = d["tissue_weight"].sum()
    rows_b = []
    for cat in GENO_ORDER:
        cat_w = d.loc[d["genotype_cat"] == cat, "tissue_weight"].sum()
        pct = (cat_w / total_w) * 100.0 if total_w > 0 else 0.0
        rows_b.append({"category": cat, "percentage": pct, "label": f"{pct:.0f}%"})
    df_b = pd.DataFrame(rows_b)
    df_b.to_csv(out_dir / "Figure4b_genotype_composition.csv", index=False)
    print("\n[Figure 4b Summary - BRCA1 Genotype Breakdown (Global)]:")
    for r in rows_b:
        print(f"  {r['category']:<25}: {r['percentage']:.2f}%")

    rows_b_age = []
    for ag in AGE_ORDER:
        sub = d[d["age_group"] == ag]
        ag_w = sub["tissue_weight"].sum()
        pm_w = sub.loc[sub["genotype_cat"] == "Exclusively BRCA1+/-", "tissue_weight"].sum()
        mix_w = sub.loc[sub["genotype_cat"] == "Mixed", "tissue_weight"].sum()
        mm_w = sub.loc[sub["genotype_cat"] == "Exclusively BRCA1-/-", "tissue_weight"].sum()
        rows_b_age.append({
            "age_group": ag,
            "exclusively_brca1_pm_pct": round((pm_w / ag_w * 100.0) if ag_w > 0 else 0.0, 2),
            "mixed_pct": round((mix_w / ag_w * 100.0) if ag_w > 0 else 0.0, 2),
            "exclusively_brca1_mm_pct": round((mm_w / ag_w * 100.0) if ag_w > 0 else 0.0, 2),
        })
    df_b_age = pd.DataFrame(rows_b_age)
    df_b_age.to_csv(out_dir / "Figure4b_genotype_composition_by_age.csv", index=False)
    print("\n[Figure 4b Summary - BRCA1 Genotype Breakdown by Age Group]:")
    print(df_b_age.to_string(index=False))

    rows_c = []
    for ag in AGE_ORDER:
        sub = d[d["age_group"] == ag]
        ag_w = sub["tissue_weight"].sum()
        row = {"age_group": ag}
        for fc in FOUNDER_ORDER:
            fc_w = sub.loc[sub["founder_category"] == fc, "tissue_weight"].sum()
            row[f"founder_{fc}_pct"] = round((fc_w / ag_w * 100.0) if ag_w > 0 else 0.0, 2)
        rows_c.append(row)
    df_c = pd.DataFrame(rows_c)
    df_c.to_csv(out_dir / "Figure4c_founder_lineages_by_age.csv", index=False)
    print("\n[Figure 4c Summary - Founder Lineage Breakdown by Age Group]:")
    print(df_c.to_string(index=False))

    # Trend stats
    age_rank = {ag: idx for idx, ag in enumerate(AGE_ORDER)}
    x_rank = d["age_group"].map(age_rank).to_numpy()
    y_founders = d["founder_count"].to_numpy()
    tau, p_val = kendalltau(x_rank, y_founders)
    print(f"\n[Statistical Validation]:")
    print(f"  Kendall's tau-b: {tau:.4f} (P = {p_val:.2e})")

    # Contingency Table
    ctable = pd.crosstab(d["age_group"], d["founder_category"]).reindex(index=AGE_ORDER, columns=FOUNDER_ORDER).fillna(0)
    ct_clean = ctable.loc[(ctable.sum(axis=1) > 0), (ctable.sum(axis=0) > 0)]
    if not ct_clean.empty and ct_clean.shape[0] > 1 and ct_clean.shape[1] > 1:
        chi2, chi2_p, dof, _ = chi2_contingency(ct_clean)
    else:
        chi2, chi2_p, dof = float("nan"), float("nan"), 0
    print(f"  Chi-square: {chi2:.2f}, df: {dof} (P = {chi2_p:.2e})")

    with open(out_dir / "stats_validation.txt", "w") as f:
        f.write(f"Kendall_tau: {tau}\nKendall_p: {p_val}\nChi2: {chi2}\nChi2_p: {chi2_p}\ndof: {dof}\n")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--post-csv", type=str, default="experiments/exp08_cs_abc_smc_homeostasis/data/generation_14.csv")
    parser.add_argument("--binary", type=str, default="./build/run_paper_features")
    parser.add_argument("--n-per-particle", type=int, default=50)
    parser.add_argument("--workers", type=int, default=10)
    parser.add_argument("--output-dir", type=str, default="experiments/exp08_cs_abc_smc_homeostasis/emergent_features")
    args = parser.parse_args()

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    raw_df = simulate_cohort(Path(args.post_csv), args.binary, args.n_per_particle, args.workers, out_dir)
    summarize_features(raw_df, out_dir)


if __name__ == "__main__":
    main()
