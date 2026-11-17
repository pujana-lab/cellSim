#!/usr/bin/env python3
"""
Unit tests for the analytical beta_prior function in abc_smc.py.
Verifies boundary conditions, monotonicity, numerical stability,
and exact satisfaction of the analytical tissue retention formula:
    lo <= ((1 - beta) * (1 - mu) * (1 + d))^T <= hi
"""

import sys
import unittest
from pathlib import Path

# Add project root and exp08 directory to sys.path
SCRIPT_DIR = Path(__file__).resolve().parent
EXP08_DIR = SCRIPT_DIR.parent / "experiments" / "exp08_cs_abc_smc_homeostasis"
sys.path.insert(0, str(EXP08_DIR))

from abc_smc import beta_prior, is_valid_particle, sample_prior, INDEPENDENT_PRIORS, PARAM_NAMES
import numpy as np


class TestBetaPrior(unittest.TestCase):

    def test_analytical_retention_bounds(self):
        """Verifies that any beta sampled within beta_prior satisfies lo <= E[B(T)]/B(0) <= hi exactly."""
        test_cases = [
            # (mu, d, lo, hi, T)
            (0.003, 0.052, 0.15, 0.35, 80),
            (0.003, 0.105, 0.15, 0.35, 80),
            (0.003, 0.140, 0.20, 0.30, 80),
            (0.001, 0.080, 0.10, 0.40, 80),
            (0.005, 0.060, 0.15, 0.35, 80),
        ]
        for mu, d, lo, hi, T in test_cases:
            b_min, b_max = beta_prior(mu, d, lo, hi, T)
            k = (1.0 - mu) * (1.0 + d)

            # At b_min, retention must equal hi
            retention_at_b_min = ((1.0 - b_min) * k) ** T
            self.assertAlmostEqual(retention_at_b_min, hi, places=6,
                                   msg=f"Failed at b_min for mu={mu}, d={d}")

            # At b_max, retention must equal lo
            retention_at_b_max = ((1.0 - b_max) * k) ** T
            self.assertAlmostEqual(retention_at_b_max, lo, places=6,
                                   msg=f"Failed at b_max for mu={mu}, d={d}")

            # Any intermediate beta must stay strictly between lo and hi
            b_mid = 0.5 * (b_min + b_max)
            retention_mid = ((1.0 - b_mid) * k) ** T
            self.assertTrue(lo <= retention_mid <= hi)

    def test_monotonicity_with_division_rate(self):
        """Higher division rates d require higher beta to maintain the same retention window."""
        mu = 0.003
        lo, hi = 0.15, 0.35
        d_values = [0.06, 0.08, 0.10, 0.12, 0.14]
        prev_b_min = -1.0
        for d in d_values:
            b_min, b_max = beta_prior(mu, d, lo, hi, 80)
            self.assertLess(b_min, b_max)
            if prev_b_min > 0:
                self.assertGreater(b_min, prev_b_min, "Beta min should strictly increase with higher d")
            prev_b_min = b_min

    def test_is_valid_particle_validation(self):
        """Checks particle validation correctly enforces conditioned beta."""
        tp53_rate = 0.003
        d_val = 0.10
        b_min, b_max = beta_prior(tp53_rate, d_val, 0.15, 0.35, 80)

        # Valid particle
        valid_theta = np.array([d_val, 0.5 * (b_min + b_max), 0.10, 0.30, 0.25])
        self.assertTrue(is_valid_particle(valid_theta, tp53_rate, 0.15, 0.35))

        # Beta too low (insufficient cell loss, tissue would over-proliferate)
        invalid_low_beta = np.array([d_val, b_min - 0.02, 0.10, 0.30, 0.25])
        self.assertFalse(is_valid_particle(invalid_low_beta, tp53_rate, 0.15, 0.35))

        # Beta too high (excessive loss, tissue decays too fast)
        invalid_high_beta = np.array([d_val, b_max + 0.02, 0.10, 0.30, 0.25])
        self.assertFalse(is_valid_particle(invalid_high_beta, tp53_rate, 0.15, 0.35))

        # Delta order violation (delta_high <= delta_low)
        invalid_delta = np.array([d_val, 0.5 * (b_min + b_max), 0.20, 0.15, 0.25])
        self.assertFalse(is_valid_particle(invalid_delta, tp53_rate, 0.15, 0.35))

    def test_sample_prior_generates_valid_particles(self):
        """Tests that 100 randomly sampled particles all satisfy the conditioned constraints."""
        rng = np.random.default_rng(123)
        tp53_rate = 0.003
        lo, hi = 0.15, 0.35
        for _ in range(100):
            p = sample_prior(rng, tp53_rate, lo, hi)
            self.assertTrue(is_valid_particle(p, tp53_rate, lo, hi))
            d = p[0]
            beta = p[1]
            b_lo, b_hi = beta_prior(tp53_rate, d, lo, hi, 80)
            self.assertTrue(b_lo <= beta <= b_hi)


if __name__ == "__main__":
    unittest.main()
