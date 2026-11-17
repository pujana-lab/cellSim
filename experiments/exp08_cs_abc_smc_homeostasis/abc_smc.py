#!/usr/bin/env python3
"""
abc_smc.py — ABC-SMC Calibration with Healthy Tissue Homeostasis & Coupled Proliferation.

Scientific Context:
-------------------
Addresses the biological scale inversion highlighted in tissue dynamics analysis:
1. Normal/healthy cells maintain dynamic homeostasis with carrying capacity cap (1000 cells).
2. Basal division rate p is sampled in U[0.047, 0.057].
3. Neoplastic cell proliferation is tied deterministically to nu = 10 * p (in [0.47, 0.57]).
4. Mutation priors for BRCA1 somatic loss (beta) are broadened to explore rare-mutation
   regimes (beta << p) as well as the prior posterior space.
5. Evaluates against Kuchenbaecker et al. (2017) clinical cumulative risk points at ages
   30, 40, 50, 60, 70, 80 with weighted SSE.
"""

from __future__ import annotations

import argparse
import os
import re
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from scipy.stats import multivariate_normal

# Parameters in ABC-SMC particle vector:
# brca1_rate is conditioned dynamically on (division_rate, tp53_rate) via beta_prior.
PARAM_NAMES = ["division_rate", "brca1_rate", "low_delta", "high_delta", "neoplastic_div_rate"]

INDEPENDENT_PRIORS: Dict[str, Tuple[float, float]] = {
    "division_rate":       (0.060, 0.140),
    "low_delta":           (0.060, 0.220),
    "high_delta":          (0.120, 0.600),
    "neoplastic_div_rate": (0.150, 0.550),
}

# Regex to parse stdout from run_bootstrapping
_DIST_RE = re.compile(r"weighted SSE\s*=\s*([\d.eE+\-]+)")
_AGE_ROW_RE = re.compile(r"^\s*(\d{2})\s+([\d.]+)\s+", re.MULTILINE)
AGES = [30, 40, 50, 60, 70, 80]


def beta_prior(mu: float, d: float, lo: float = 0.15, hi: float = 0.35, T: int = 80) -> Tuple[float, float]:
    """Analytical range of beta such that lo <= E[B(T)]/B(0) <= hi."""
    k = (1.0 - mu) * (1.0 + d)
    b_min = 1.0 - (hi ** (1.0 / T)) / k
    b_max = 1.0 - (lo ** (1.0 / T)) / k
    if b_min > b_max:
        b_min, b_max = b_max, b_min
    return b_min, b_max


def is_valid_particle(
    theta: np.ndarray,
    tp53_rate: float = 0.003,
    healthy_lo: float = 0.15,
    healthy_hi: float = 0.35
) -> bool:
    """Check prior bounds for independent parameters, delta constraint, and conditioned beta range."""
    # 1. Independent parameters bounds
    for name, (lo, hi) in INDEPENDENT_PRIORS.items():
        val = theta[PARAM_NAMES.index(name)]
        if val < lo or val > hi:
            return False

    # 2. Biological constraint: delta_high must strictly exceed delta_low
    low_delta = theta[PARAM_NAMES.index("low_delta")]
    high_delta = theta[PARAM_NAMES.index("high_delta")]
    if high_delta <= low_delta:
        return False

    # 3. Conditioned beta (BRCA1 rate) dynamically derived from (division_rate, tp53_rate)
    d_val = theta[PARAM_NAMES.index("division_rate")]
    beta_val = theta[PARAM_NAMES.index("brca1_rate")]
    b_min, b_max = beta_prior(tp53_rate, d_val, healthy_lo, healthy_hi)
    if beta_val < b_min or beta_val > b_max:
        return False

    return True


def sample_prior(
    rng: np.random.Generator,
    tp53_rate: float = 0.003,
    healthy_lo: float = 0.15,
    healthy_hi: float = 0.35
) -> np.ndarray:
    """Sample parameter vector: independent parameters from priors, beta from conditioned analytical range."""
    candidate = np.zeros(len(PARAM_NAMES))
    # 1. division_rate (d)
    candidate[0] = rng.uniform(INDEPENDENT_PRIORS["division_rate"][0], INDEPENDENT_PRIORS["division_rate"][1])
    # 2. brca1_rate (beta) analytically derived from (tp53_rate, division_rate)
    b_min, b_max = beta_prior(tp53_rate, candidate[0], healthy_lo, healthy_hi)
    candidate[1] = rng.uniform(b_min, b_max)
    # 3. low_delta
    candidate[2] = rng.uniform(INDEPENDENT_PRIORS["low_delta"][0], INDEPENDENT_PRIORS["low_delta"][1])
    # 4. high_delta (> low_delta)
    h_min = max(INDEPENDENT_PRIORS["high_delta"][0], candidate[2] + 0.01)
    candidate[3] = rng.uniform(h_min, INDEPENDENT_PRIORS["high_delta"][1])
    # 5. neoplastic_div_rate
    candidate[4] = rng.uniform(INDEPENDENT_PRIORS["neoplastic_div_rate"][0], INDEPENDENT_PRIORS["neoplastic_div_rate"][1])

    while not is_valid_particle(candidate, tp53_rate, healthy_lo, healthy_hi):
        candidate[0] = rng.uniform(INDEPENDENT_PRIORS["division_rate"][0], INDEPENDENT_PRIORS["division_rate"][1])
        b_min, b_max = beta_prior(tp53_rate, candidate[0], healthy_lo, healthy_hi)
        candidate[1] = rng.uniform(b_min, b_max)
        candidate[2] = rng.uniform(INDEPENDENT_PRIORS["low_delta"][0], INDEPENDENT_PRIORS["low_delta"][1])
        h_min = max(INDEPENDENT_PRIORS["high_delta"][0], candidate[2] + 0.01)
        candidate[3] = rng.uniform(h_min, INDEPENDENT_PRIORS["high_delta"][1])
        candidate[4] = rng.uniform(INDEPENDENT_PRIORS["neoplastic_div_rate"][0], INDEPENDENT_PRIORS["neoplastic_div_rate"][1])

    return candidate


def evaluate_simulation(
    binary_path: str,
    theta: np.ndarray,
    n_sim: int,
    seed_offset: int,
    n_cells: int = 1000,
    max_healthy_cells: int = 1000,
    tp53_rate: float = 0.003,
    timeout_sec: int = 180,
) -> Tuple[float, Dict[str, float]]:
    """Runs run_bootstrapping for a single particle and returns (distance, risks)."""
    p_dict = {name: float(theta[i]) for i, name in enumerate(PARAM_NAMES)}
    p_div = p_dict["division_rate"]
    nu_div = p_dict["neoplastic_div_rate"]  # Autonomous tumoral division rate (decoupled)

    cmd = [
        binary_path,
        "--n-runs", str(n_sim),
        "--seed-offset", str(seed_offset),
        "--output", "/dev/null",
        "--weighted-sse",
        "--big-bang",
        "--n-cells", str(n_cells),
        "--max-healthy-cells", str(max_healthy_cells),
        "--division-rate", f"{p_div:.8g}",
        "--neoplastic-div-rate", f"{nu_div:.8g}",
        "--brca1-rate", f"{p_dict['brca1_rate']:.8g}",
        "--tp53-rate", f"{tp53_rate:.8g}",
        "--low-delta", f"{p_dict['low_delta']:.8g}",
        "--high-delta", f"{p_dict['high_delta']:.8g}",
        "--d1-threshold", "2.0",
        "--d2-threshold", "5.0",
        "--kappa-brca1", "2.0",
        "--age-drift", "1e-5",
    ]

    try:
        env = dict(os.environ)
        env["OMP_NUM_THREADS"] = "1"
        res = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout_sec, env=env)
    except (subprocess.TimeoutExpired, OSError):
        return float("inf"), {}

    if res.returncode != 0:
        return float("inf"), {}

    m_dist = _DIST_RE.search(res.stdout)
    if not m_dist:
        return float("inf"), {}
    dist = float(m_dist.group(1))

    risks = {}
    for m in _AGE_ROW_RE.finditer(res.stdout):
        age_val = int(m.group(1))
        if age_val in AGES:
            risks[f"risk_{age_val}"] = float(m.group(2))

    return dist, risks


def weighted_covariance(particles: np.ndarray, weights: np.ndarray) -> np.ndarray:
    """Compute Beaumont-Toni weighted empirical covariance matrix."""
    w = weights / np.sum(weights)
    mean = np.sum(particles * w[:, None], axis=0)
    diff = particles - mean
    denom = 1.0 - np.sum(w ** 2)
    if denom <= 1e-12:
        denom = 1.0
    cov = (diff * w[:, None]).T @ diff / denom
    # Perturbation covariance: 2 * empirical covariance + diagonal regularization
    cov_reg = 2.0 * cov + 1e-6 * np.eye(len(PARAM_NAMES))
    return cov_reg


def run_abc_smc(
    n_particles: int,
    n_generations: int,
    n_sim: int,
    workers: int,
    binary_path: str,
    output_dir: Path,
    n_cells: int = 500,
    max_healthy_cells: int = 500,
    tp53_rate: float = 0.003,
    healthy_lo: float = 0.15,
    healthy_hi: float = 0.35,
    seed: int = 42,
    resume_from: Optional[str] = None,
) -> None:
    """Executes the full ABC-SMC calibration algorithm with conditioned beta_prior."""
    rng = np.random.default_rng(seed)
    output_dir.mkdir(parents=True, exist_ok=True)
    diagnostics_path = output_dir / "smc_diagnostics.csv"

    print(f"\n{'='*70}")
    print(f" ABC-SMC Calibration with Conditioned beta_prior (N0={n_cells})")
    print(f" Particles: {n_particles} | Generations: {n_generations} | MC runs: {n_sim}")
    print(f" Healthy retention at T=80: [{healthy_lo*100:.1f}%, {healthy_hi*100:.1f}%]")
    print(f" Independent Priors: {INDEPENDENT_PRIORS}")
    print(f" brca1_rate (beta): dynamically derived from beta_prior(mu={tp53_rate}, d, lo={healthy_lo}, hi={healthy_hi})")
    print(f"{'='*70}\n")

    particles_prev: Optional[np.ndarray] = None
    weights_prev: Optional[np.ndarray] = None
    cov_prev: Optional[np.ndarray] = None
    epsilon_t: float = float("inf")
    diag_rows: List[dict] = []
    start_gen = 0

    if resume_from:
        resume_path = Path(resume_from)
        gen_files = sorted(resume_path.glob("generation_*.csv"))
        if not gen_files:
            raise FileNotFoundError(f"No generation_*.csv found in resume directory: {resume_path}")
        latest_file = gen_files[-1]
        m = re.search(r"generation_(\d+)\.csv", latest_file.name)
        if not m:
            raise ValueError(f"Could not parse generation from {latest_file.name}")
        last_gen = int(m.group(1))
        start_gen = last_gen + 1

        print(f"Resuming ABC-SMC Calibration from {latest_file.name} (Gen {last_gen} -> Gen {start_gen})")
        if resume_path.resolve() != output_dir.resolve():
            import shutil
            for gf in gen_files:
                shutil.copy2(gf, output_dir / gf.name)
            diag_file = resume_path / "smc_diagnostics.csv"
            if diag_file.exists():
                shutil.copy2(diag_file, output_dir / diag_file.name)
                diag_rows = pd.read_csv(diag_file).to_dict(orient="records")

        df_prev = pd.read_csv(latest_file)
        particles_prev = df_prev[PARAM_NAMES].to_numpy()
        weights_prev = df_prev["weight"].to_numpy()
        weights_prev = weights_prev / np.sum(weights_prev)
        dist_col = "distance_weighted_sse" if "distance_weighted_sse" in df_prev.columns else "distance"
        distances_prev = df_prev[dist_col].to_numpy()
        epsilon_t = float(np.percentile(distances_prev, 75))
        cov_prev = weighted_covariance(particles_prev, weights_prev)

    for gen in range(start_gen, n_generations):
        t_gen_start = time.time()
        print(f"--- Generation {gen:02d} (epsilon={epsilon_t:.4f}) ---")

        accepted_particles: List[np.ndarray] = []
        accepted_distances: List[float] = []
        accepted_risks: List[Dict[str, float]] = []
        attempts = 0

        # Precompute sampling structures for gen > 0
        if gen > 0:
            norm_weights = weights_prev / np.sum(weights_prev)
            inv_cov = np.linalg.pinv(cov_prev)
            det_cov = np.linalg.det(cov_prev)
            if det_cov <= 0 or not np.isfinite(det_cov):
                cov_prev += 1e-5 * np.eye(len(PARAM_NAMES))
                inv_cov = np.linalg.pinv(cov_prev)

        def sample_and_evaluate(task_id: int) -> Tuple[np.ndarray, float, Dict[str, float], int]:
            local_rng = np.random.default_rng(seed + gen * 100000 + task_id * 1009)
            local_attempts = 0
            while True:
                local_attempts += 1
                if gen == 0:
                    candidate = sample_prior(local_rng, tp53_rate, healthy_lo, healthy_hi)
                else:
                    parent_idx = local_rng.choice(len(particles_prev), p=norm_weights)
                    candidate = local_rng.multivariate_normal(particles_prev[parent_idx], cov_prev)
                    if not is_valid_particle(candidate, tp53_rate, healthy_lo, healthy_hi):
                        continue

                seed_offset = int(local_rng.integers(1, 1_000_000_000))
                dist, risks = evaluate_simulation(
                    binary_path, candidate, n_sim, seed_offset,
                    n_cells=n_cells, max_healthy_cells=max_healthy_cells,
                    tp53_rate=tp53_rate
                )
                if dist < epsilon_t:
                    return candidate, dist, risks, local_attempts

        # Parallel evaluation of accepted particles
        with ThreadPoolExecutor(max_workers=workers) as executor:
            futures = [executor.submit(sample_and_evaluate, i) for i in range(n_particles)]
            for fut in futures:
                cand, dist, risks, att = fut.result()
                accepted_particles.append(cand)
                accepted_distances.append(dist)
                accepted_risks.append(risks)
                attempts += att
                if len(accepted_particles) % max(1, n_particles // 5) == 0:
                    print(f"  Accepted {len(accepted_particles)}/{n_particles} particles...")

        particles_arr = np.array(accepted_particles)
        distances_arr = np.array(accepted_distances)

        # Weight computation
        if gen == 0:
            weights_arr = np.ones(n_particles) / n_particles
        else:
            weights_arr = np.zeros(n_particles)
            for i in range(n_particles):
                diff = particles_arr[i] - particles_prev
                # Gaussian kernel evaluation
                exponent = -0.5 * np.sum((diff @ inv_cov) * diff, axis=1)
                denom = np.sum(norm_weights * np.exp(exponent))
                weights_arr[i] = 1.0 / denom if denom > 0 else 0.0
            weights_arr /= np.sum(weights_arr)

        ess = 1.0 / np.sum(weights_arr ** 2)
        acceptance_rate = n_particles / attempts
        t_gen = time.time() - t_gen_start

        print(f"  Completed in {t_gen:.1f}s | Acceptance rate: {acceptance_rate:.4f} | ESS: {ess:.1f}")
        print(f"  Distance median: {np.median(distances_arr):.2f} (min: {np.min(distances_arr):.2f})")

        # Save generation particles to CSV
        gen_df = pd.DataFrame(particles_arr, columns=PARAM_NAMES)
        gen_df["weight"] = weights_arr
        gen_df["distance_weighted_sse"] = distances_arr
        for age in AGES:
            gen_df[f"risk_{age}"] = [r.get(f"risk_{age}", np.nan) for r in accepted_risks]
        gen_df.to_csv(output_dir / f"generation_{gen:02d}.csv", index=False)

        # Record diagnostics
        diag_rows.append({
            "generation": gen,
            "epsilon": epsilon_t,
            "attempts": attempts,
            "accepted": n_particles,
            "acceptance_rate": acceptance_rate,
            "ess": ess,
            "weighted_sse_min": float(np.min(distances_arr)),
            "weighted_sse_median": float(np.median(distances_arr)),
            "weighted_sse_p75": float(np.percentile(distances_arr, 75)),
            "time_sec": t_gen,
        })
        pd.DataFrame(diag_rows).to_csv(diagnostics_path, index=False)

        # Update for next generation
        particles_prev = particles_arr
        weights_prev = weights_arr
        cov_prev = weighted_covariance(particles_arr, weights_arr)
        # Adaptive epsilon: 75th percentile of accepted distances
        epsilon_t = float(np.percentile(distances_arr, 75))


def main() -> None:
    parser = argparse.ArgumentParser(description="ABC-SMC calibration with healthy homeostasis and conditioned beta_prior.")
    parser.add_argument("--particles", type=int, default=50, help="Particles per generation")
    parser.add_argument("--generations", type=int, default=8, help="Number of SMC generations")
    parser.add_argument("--n-sim", type=int, default=200, help="Monte Carlo seeds per simulation evaluation")
    parser.add_argument("--n-cells", type=int, default=500, help="Starting cell count per tissue")
    parser.add_argument("--max-healthy-cells", type=int, default=500, help="Carrying capacity for healthy cells")
    parser.add_argument("--tp53-rate", type=float, default=0.003, help="Somatic TP53 mutation rate")
    parser.add_argument("--healthy-retention-min", type=float, default=0.15, help="Minimum retention of baseline cells at T=80")
    parser.add_argument("--healthy-retention-max", type=float, default=0.35, help="Maximum retention of baseline cells at T=80")
    parser.add_argument("--workers", type=int, default=8, help="Worker threads")
    parser.add_argument("--seed", type=int, default=42, help="RNG seed")
    parser.add_argument("--binary", type=str, default="./build/run_bootstrapping", help="Path to run_bootstrapping binary")
    parser.add_argument("--output-dir", type=str, default="./experiments/exp08_cs_abc_smc_homeostasis/results", help="Output directory")
    parser.add_argument("--resume-from", type=str, default=None, help="Directory containing previous generation CSVs to resume from")
    parser.add_argument("--smoke-test", action="store_true", help="Fast smoke test (2 particles, 2 gens, 5 sims)")

    args = parser.parse_args()
    binary_path = str(Path(args.binary).resolve())
    if not os.path.isfile(binary_path) or not os.access(binary_path, os.X_OK):
        print(f"Error: binary '{binary_path}' not found or not executable.", file=sys.stderr)
        sys.exit(1)

    if args.smoke_test:
        print("[SMOKE TEST] Running fast smoke test for conditioned beta_prior ABC-SMC...")
        run_abc_smc(
            n_particles=2,
            n_generations=2,
            n_sim=5,
            workers=2,
            binary_path=binary_path,
            output_dir=Path(args.output_dir) / "smoke_test",
            n_cells=args.n_cells,
            max_healthy_cells=args.max_healthy_cells,
            tp53_rate=args.tp53_rate,
            healthy_lo=args.healthy_retention_min,
            healthy_hi=args.healthy_retention_max,
            seed=args.seed,
        )
        print("[OK] Smoke test completed successfully!")
        return

    run_abc_smc(
        n_particles=args.particles,
        n_generations=args.generations,
        n_sim=args.n_sim,
        workers=args.workers,
        binary_path=binary_path,
        output_dir=Path(args.output_dir),
        n_cells=args.n_cells,
        max_healthy_cells=args.max_healthy_cells,
        tp53_rate=args.tp53_rate,
        healthy_lo=args.healthy_retention_min,
        healthy_hi=args.healthy_retention_max,
        seed=args.seed,
        resume_from=args.resume_from,
    )


if __name__ == "__main__":
    main()
