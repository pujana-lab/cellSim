#!/usr/bin/env python3
"""
generate_figures_1000cells_final.py — Publication figures for 1000-cell Gen 10 Calibration:
  - Figure 4 (Posteriors, Violin plot of BRCA1-/- clonal purity, Founder complexity)
  - Figure 3 (ABC-SMC convergence to SSE 1.08, Efficiency/ESS, Penetrance vs Kuchenbaecker, Neutral dynamics)
  - Figure_stage_dynamics (Cell stage evolution, relative composition, allelic breakdown, milestone table)
"""

from __future__ import annotations

from pathlib import Path
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import kendalltau

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

SCRIPT_DIR = Path(__file__).resolve().parent
RESULTS_DIR = SCRIPT_DIR / "results_1000cells_final"
RUN_DIR = SCRIPT_DIR / "runs" / "run_abc_xi_wide_20261120_171832"

POST_CSV = RESULTS_DIR / "posterior_particles.csv"
STAGE_CSV = RESULTS_DIR / "stage_trajectories_1000.csv"
FEATURES_CSV = RESULTS_DIR / "features_1000.csv"
INTERCEPTION_DESIGN_CSV = RESULTS_DIR / "interception" / "interception_design.csv"
INTERCEPTION_RESULTS_CSV = RESULTS_DIR / "interception" / "interception_results.csv"
DIAGNOSTICS_CSV = RESULTS_DIR / "abc_diagnostics.csv"

AGE_ORDER = ["<30", "30-39", "40-49", "50-59", "60-69", "70-80"]
AGES_CLINICAL = np.array([30, 40, 50, 60, 70, 80], dtype=int)
OBS_CLINICAL = np.array([4.0, 24.0, 43.0, 56.0, 65.0, 72.0])
CI_LO = np.array([2.0, 20.0, 38.0, 51.0, 60.0, 65.0])
CI_HI = np.array([8.0, 29.0, 48.0, 61.0, 70.0, 79.0])


def age_group(onset: int) -> str:
    if onset < 30: return "<30"
    if onset < 40: return "30-39"
    if onset < 50: return "40-49"
    if onset < 60: return "50-59"
    if onset < 70: return "60-69"
    return "70-80"


def render_figure4():
    df_feat = pd.read_csv(FEATURES_CSV)
    df_onset = df_feat[df_feat["onset_year"] > 0].copy()
    df_onset["age_cat"] = df_onset["onset_year"].apply(age_group)
    df_onset["frac_mm"] = (df_onset["brca1_mm"] / df_onset["tumoral_cells"]) * 100.0

    # Panel c founder counts (percentages)
    def map_founder(f: int) -> str:
        if f == 1: return "1"
        elif f == 2: return "2"
        elif f == 3: return "3"
        else: return ">=4"

    df_onset["founder_cat"] = df_onset["founder_count"].apply(map_founder)
    ct_found_raw = pd.crosstab(df_onset["age_cat"], df_onset["founder_cat"])
    for ag in AGE_ORDER:
        if ag not in ct_found_raw.index:
            ct_found_raw.loc[ag] = 0
    ct_found_raw = ct_found_raw.reindex(AGE_ORDER).fillna(0)
    for f_col in ["1", "2", "3", ">=4"]:
        if f_col not in ct_found_raw.columns:
            ct_found_raw[f_col] = 0

    # Relative percentage table per age group
    row_sums = ct_found_raw.sum(axis=1)
    ct_found = ct_found_raw.div(row_sums.replace(0, 1), axis=0) * 100.0

    tau, p_val = kendalltau(df_onset["onset_year"], df_onset["founder_count"])

    # Load Posterior Particles for Figure 4a
    df_part = pd.read_csv(POST_CSV)
    df_part["high_delta_over_xi"] = df_part["high_delta"] / df_part["xi"]
    w = df_part["weight"].to_numpy(float)
    w = w / np.sum(w)

    fig = plt.figure(figsize=(12.0, 11.0), dpi=300)
    fig.text(0.04, 0.98, "Figure 4", fontsize=16, fontweight="bold", va="top")

    # Panel a: Posteriors
    fig.text(0.04, 0.94, "a", fontsize=14, fontweight="bold", va="top")
    fig.text(0.07, 0.94, "Calibrated parameter posteriors (1000 cells neutral dynamics, F = 1.0)", fontsize=11, fontweight="bold", va="top")

    gs_a = fig.add_gridspec(
        nrows=2, ncols=5,
        left=0.07, right=0.96, top=0.90, bottom=0.56,
        width_ratios=[1.0, 1.0, 1.0, 1.0, 1.15],
        wspace=0.38, hspace=0.35
    )

    params = ["high_delta_over_xi", "low_delta", "high_delta", "neoplastic_div_rate"]
    param_labels = {
        "high_delta_over_xi": r"$\delta_{\rm high} / \xi$",
        "low_delta": r"$\delta_{\rm low}$",
        "high_delta": r"$\delta_{\rm high}$",
        "neoplastic_div_rate": r"$\nu$",
    }
    priors_range = {
        "high_delta_over_xi": (0.0, 0.035),
        "low_delta": (0.05, 0.25),
        "high_delta": (0.10, 0.50),
        "neoplastic_div_rate": (0.05, 0.60),
    }

    # Precompute prior distribution for ratio delta_high / xi via Monte Carlo
    rng_prior = np.random.default_rng(20261125)
    xi_mc = rng_prior.uniform(1.0, 50.0, 200_000)
    hd_mc = rng_prior.uniform(0.10, 0.50, 200_000)
    ratio_mc = hd_mc / xi_mc
    kde_ratio_prior = stats.gaussian_kde(ratio_mc)

    for col_idx, p_name in enumerate(params):
        ax_top = fig.add_subplot(gs_a[0, col_idx])
        ax_bot = fig.add_subplot(gs_a[1, col_idx])

        vals = df_part[p_name].to_numpy()
        lo, hi = priors_range[p_name]

        idx_sort = np.argsort(vals)
        cum_w = np.cumsum(w[idx_sort])
        q50 = float(np.interp(0.50, cum_w, vals[idx_sort]))

        xs = np.linspace(lo, hi, 200)
        kde = stats.gaussian_kde(vals, weights=w)
        ax_top.plot(xs, kde(xs), color="black", lw=1.5, label="Posterior KDE")

        if p_name == "high_delta_over_xi":
            ax_top.plot(xs, kde_ratio_prior(xs), color="#9CA3AF", linestyle="--", lw=1.0, label="Prior")
        else:
            ax_top.axhline(1.0 / (hi - lo), color="#9CA3AF", linestyle="--", lw=1.0, label="Uniform prior")

        ax_top.axvline(q50, color="#DC2626", lw=1.5, label="Median")
        fmt_q50 = f" {q50:.4f}" if p_name == "high_delta_over_xi" else f" {q50:.3f}"
        ax_top.text(q50, ax_top.get_ylim()[1] * 0.9, fmt_q50, color="#DC2626", fontweight="bold", fontsize=8.5, ha="center")
        ax_top.set_title(param_labels[p_name], fontsize=11, fontweight="bold")
        ax_top.set_xlim(lo, hi)
        if col_idx == 0:
            ax_top.set_ylabel("Density", fontsize=9.5)
        if col_idx == 3:
            ax_top.legend(fontsize=7.0, frameon=False, loc="upper right")

        ax_bot.plot(vals[idx_sort], cum_w, color="black", lw=1.5, label="Posterior CDF")
        if p_name == "high_delta_over_xi":
            prior_cdf = np.array([(ratio_mc <= x).mean() for x in xs])
            ax_bot.plot(xs, prior_cdf, color="#9CA3AF", linestyle="--", lw=1.0, label="Theoretical prior")
        else:
            ax_bot.plot([lo, hi], [0.0, 1.0], color="#9CA3AF", linestyle="--", lw=1.0, label="Theoretical prior")
        ax_bot.set_xlim(lo, hi)
        ax_bot.set_ylim(0, 1.02)
        if col_idx == 0:
            ax_bot.set_ylabel("CDF", fontsize=9.5)
        if col_idx == 3:
            ax_bot.legend(fontsize=7.0, frameon=False, loc="center right")

    # Correlation matrix
    ax_corr = fig.add_subplot(gs_a[:, 4])
    X = df_part[params].to_numpy()
    cov_X = np.cov(X, rowvar=False, aweights=w)
    std_X = np.sqrt(np.diag(cov_X))
    corr_X = cov_X / np.outer(std_X, std_X)

    im_c = ax_corr.imshow(corr_X, cmap="coolwarm", vmin=-1, vmax=1)
    ax_corr.set_xticks(range(4))
    ax_corr.set_yticks(range(4))
    ax_corr.set_xticklabels([param_labels[p] for p in params], fontsize=8.5, rotation=45, ha="right")
    ax_corr.set_yticklabels([param_labels[p] for p in params], fontsize=8.5)
    ax_corr.set_title("Weighted\ncorrelation", fontsize=9.5, fontweight="bold", pad=8)
    cb = fig.colorbar(im_c, ax=ax_corr, fraction=0.046, pad=0.06)
    cb.set_ticks([-1, 0, 1])

    # Panel b: Barplot of BRCA1 status (Pure BRCA1-/- vs Mixed BRCA1)
    fig.text(0.04, 0.50, "b", fontsize=14, fontweight="bold", va="top")
    ax_b = fig.add_axes([0.07, 0.08, 0.40, 0.38])
    ax_b.set_title("BRCA1-/- clonal status by detection age", fontsize=11, fontweight="bold", loc="left", pad=10)

    # Classify tumors: Pure BRCA1-/- (>=80% -/-) vs Mixed BRCA1 (<80% -/- or >=20% +/-)
    df_onset["pure_brca1"] = df_onset["frac_mm"] >= 80.0
    df_onset["status_cat"] = np.where(df_onset["pure_brca1"], "Pure BRCA1-/- (≥80%)", "Mixed BRCA1 (<80%)")

    ct_status_raw = pd.crosstab(df_onset["age_cat"], df_onset["status_cat"])
    for ag in AGE_ORDER:
        if ag not in ct_status_raw.index:
            ct_status_raw.loc[ag] = 0
    ct_status_raw = ct_status_raw.reindex(AGE_ORDER).fillna(0)
    for s_col in ["Pure BRCA1-/- (≥80%)", "Mixed BRCA1 (<80%)"]:
        if s_col not in ct_status_raw.columns:
            ct_status_raw[s_col] = 0

    # Relative percentage per age group
    row_sums_b = ct_status_raw.sum(axis=1)
    ct_status = ct_status_raw.div(row_sums_b.replace(0, 1), axis=0) * 100.0

    x_positions = np.arange(len(AGE_ORDER))
    pure_vals = ct_status["Pure BRCA1-/- (≥80%)"].values
    mixed_vals = ct_status["Mixed BRCA1 (<80%)"].values

    bar_pure = ax_b.bar(
        x_positions, pure_vals, width=0.55,
        label="Pure BRCA1-/- (≥80%)", color="#2563EB", edgecolor="#111827", linewidth=0.5
    )
    bar_mixed = ax_b.bar(
        x_positions, mixed_vals, width=0.55, bottom=pure_vals,
        label="Mixed BRCA1 (<80%)", color="#93C5FD", edgecolor="#111827", linewidth=0.5
    )

    ax_b.set_xticks(x_positions)
    ax_b.set_xticklabels(AGE_ORDER, fontsize=9)
    ax_b.set_xlabel("Age at first tumor-positive detection (years)", fontsize=10)
    ax_b.set_ylabel("Proportion of tumors (%)", fontsize=10)
    ax_b.set_ylim(0, 108)
    ax_b.grid(True, linestyle="--", alpha=0.3, axis="y")
    ax_b.legend(title="Tumoral BRCA1 status", fontsize=8.0, title_fontsize=8.5, frameon=True, loc="lower right", edgecolor="#E5E7EB")

    # Panel c: Founder complexity (Relative percentage)
    fig.text(0.52, 0.50, "c", fontsize=14, fontweight="bold", va="top")
    ax_c = fig.add_axes([0.55, 0.08, 0.40, 0.38])
    ax_c.set_title("Founder-lineage complexity by detection age", fontsize=11, fontweight="bold", loc="left", pad=10)

    founder_colors = {"1": "#F59E0B", "2": "#EA580C", "3": "#DC2626", ">=4": "#7F1D1D"}
    bottom_f = np.zeros(len(AGE_ORDER))
    for f_cat in ["1", "2", "3", ">=4"]:
        vals_f = ct_found[f_cat].values
        ax_c.bar(x_positions, vals_f, width=0.55, bottom=bottom_f, label=f_cat, color=founder_colors[f_cat], edgecolor="#111827", linewidth=0.5)
        bottom_f += vals_f

    ax_c.set_xticks(x_positions)
    ax_c.set_xticklabels(AGE_ORDER, fontsize=9)
    ax_c.set_xlabel("Age at first tumor-positive detection (years)", fontsize=10)
    ax_c.set_ylabel("Proportion of tumors (%)", fontsize=10)
    ax_c.set_ylim(0, 108)
    ax_c.grid(True, linestyle="--", alpha=0.3, axis="y")
    ax_c.legend(title="Tumoral founder lineages (n)", fontsize=8.0, title_fontsize=8.5, frameon=True, loc="upper right", edgecolor="#E5E7EB")

    ax_c.text(0.50, 0.93, f"Kendall's $\\tau_b = {tau:.2f}$ ($P = {p_val:.2e}$)",
              transform=ax_c.transAxes, fontsize=8.5, fontstyle="italic", ha="center")

    out_fig4_png = RESULTS_DIR / "Figure4.png"
    out_fig4_pdf = RESULTS_DIR / "Figure4.pdf"
    plt.savefig(out_fig4_png, dpi=300, bbox_inches="tight")
    plt.savefig(out_fig4_pdf, bbox_inches="tight")
    plt.close()
    print(f"Figure 4 guardada en:\n  -> {out_fig4_png}")


def render_figure3():
    # Diagnostics across generations from cached CSV or computed
    if DIAGNOSTICS_CSV.exists():
        df_diag = pd.read_csv(DIAGNOSTICS_CSV)
    else:
        tot_eval = [200, 200, 294, 364, 726, 668, 1364, 2378, 4716, 12818]
        diag_rows = []
        for g in range(1, 11):
            f = RUN_DIR / f"particles_gen_{g}.csv"
            df_g = pd.read_csv(f)
            sse = df_g["weighted_sse"].values
            w = df_g["weight"].values
            w_norm = w / np.sum(w)
            ess = 1.0 / np.sum(w_norm ** 2)
            acc_rate = 100.0 / tot_eval[g - 1]

            diag_rows.append({
                "gen": g,
                "min_sse": np.min(sse),
                "q25_sse": np.percentile(sse, 25),
                "med_sse": np.median(sse),
                "q75_sse": np.percentile(sse, 75),
                "acc_rate": acc_rate * 100.0,
                "ess_pct": (ess / 100.0) * 100.0,
            })
        df_diag = pd.DataFrame(diag_rows)

    df_stage = pd.read_csv(STAGE_CSV)

    fig = plt.figure(figsize=(12.0, 9.8), dpi=300)
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.05], hspace=0.35, wspace=0.28)
    fig.text(0.04, 0.98, "Figure 3", fontsize=15, fontweight="bold", va="top")

    # Panel a: Weighted SSE
    ax_a = fig.add_subplot(gs[0, 0])
    ax_a.text(-0.16, 1.12, "a", transform=ax_a.transAxes, fontsize=15, fontweight="bold", va="top")
    ax_a.set_title("ABC-SMC calibration performance", fontsize=11, fontweight="bold", pad=12)

    gens = df_diag["gen"].values
    ax_a.fill_between(gens, df_diag["min_sse"], df_diag["q75_sse"], color="#E5E7EB", alpha=0.8, label="IQR / Min range")
    ax_a.plot(gens, df_diag["med_sse"], color="#DC2626", marker="o", markersize=4.5, lw=1.6, label="Median weighted SSE")
    ax_a.plot(gens, df_diag["min_sse"], color="#16A34A", linestyle="--", marker="s", markersize=4.0, lw=1.3, label="Minimum weighted SSE")
    ax_a.set_yscale("log")
    ax_a.set_xlabel("SMC Generation", fontsize=9.5)
    ax_a.set_ylabel("Weighted SSE (log scale)", fontsize=9.5)
    ax_a.set_xticks(gens)
    ax_a.grid(True, linestyle="--", alpha=0.4)
    ax_a.legend(loc="upper right", fontsize=8.0, frameon=True, facecolor="white", edgecolor="#D1D5DB")

    # Panel b: Sampling efficiency & ESS
    ax_b = fig.add_subplot(gs[0, 1])
    ax_b.text(-0.16, 1.12, "b", transform=ax_b.transAxes, fontsize=15, fontweight="bold", va="top")
    ax_b.set_title("Sampling efficiency and ESS", fontsize=11, fontweight="bold", pad=12)
    ax_b.plot(gens, df_diag["acc_rate"], color="#2563EB", marker="^", markersize=5.0, lw=1.5, label="Acceptance rate (%)")
    ax_b.plot(gens, df_diag["ess_pct"], color="#4B5563", linestyle="-.", marker="d", markersize=4.5, lw=1.4, label="ESS (% of particles)")
    ax_b.set_xlabel("SMC Generation", fontsize=9.5)
    ax_b.set_ylabel("Efficiency metric (%)", fontsize=9.5)
    ax_b.set_xticks(gens)
    ax_b.set_ylim(-2, 105)
    ax_b.grid(True, linestyle="--", alpha=0.4)
    ax_b.legend(loc="upper right", fontsize=8.0, frameon=True, facecolor="white", edgecolor="#D1D5DB")

    # Panel c: Cumulative Penetrance (Median posterior + 95% CI vs Kuchenbaecker)
    ax_c = fig.add_subplot(gs[1, 0])
    ax_c.text(-0.16, 1.12, "c", transform=ax_c.transAxes, fontsize=15, fontweight="bold", va="top")
    ax_c.set_title("Cumulative penetrance vs. clinical benchmark", fontsize=11, fontweight="bold", pad=12)

    # Compute penetrance curves per posterior particle from control tissues (condition == 'none')
    df_des = pd.read_csv(INTERCEPTION_DESIGN_CSV)
    df_res = pd.read_csv(INTERCEPTION_RESULTS_CSV)
    none_des = df_des[df_des["condition"] == "none"].copy()
    none_res = df_res[df_res["condition"] == "none"].copy()
    merged = pd.merge(none_des, none_res[["tissue_id", "onset_year"]], on="tissue_id")

    # 1. Parameter uncertainty: Posterior Credible Interval across particles
    param_cols = ["brca1_rate", "tp53_rate", "division_rate", "low_delta", "high_delta", "neoplastic_div_rate"]
    grouped = merged.groupby(param_cols)
    particle_curves = []
    for _, grp in grouped:
        c = [(grp["onset_year"].between(1, a)).mean() * 100 for a in AGES_CLINICAL]
        particle_curves.append(c)
    particle_curves = np.array(particle_curves)

    med_penetrance = np.median(particle_curves, axis=0)
    cri_lo_penetrance = np.percentile(particle_curves, 2.5, axis=0)
    cri_hi_penetrance = np.percentile(particle_curves, 97.5, axis=0)

    # 2. Simulation sampling uncertainty: Bootstrap 95% CI across the N=3000 control tissues
    rng = np.random.default_rng(20261125)
    n_tissues = len(merged)
    boot_curves = []
    for _ in range(1000):
        sample_onsets = merged["onset_year"].sample(n=n_tissues, replace=True, random_state=rng)
        boot_curves.append([(sample_onsets.between(1, a)).mean() * 100 for a in AGES_CLINICAL])
    boot_curves = np.array(boot_curves)
    boot_ci_lo = np.percentile(boot_curves, 2.5, axis=0)
    boot_ci_hi = np.percentile(boot_curves, 97.5, axis=0)

    # Shaded Posterior 95% CrI (parameter uncertainty)
    ax_c.fill_between(
        AGES_CLINICAL, cri_lo_penetrance, cri_hi_penetrance,
        color="#F87171", alpha=0.25, label="Posterior 95% CrI (params)"
    )

    # Shaded Simulation 95% CI (bootstrap across 3,000 cases)
    ax_c.fill_between(
        AGES_CLINICAL, boot_ci_lo, boot_ci_hi,
        color="#DC2626", alpha=0.45, label="Simulation 95% CI (bootstrap N=3,000)"
    )

    # Posterior Median curve
    ax_c.plot(
        AGES_CLINICAL, med_penetrance,
        color="#B91C1C", lw=2.2, marker="o", markersize=4.5,
        label="Posterior median"
    )

    # Kuchenbaecker et al. (2017) clinical benchmark
    yerr = np.vstack([OBS_CLINICAL - CI_LO, CI_HI - OBS_CLINICAL])
    ax_c.errorbar(
        AGES_CLINICAL, OBS_CLINICAL, yerr=yerr, fmt="s", color="#111827", ecolor="#111827",
        elinewidth=1.3, capsize=4.5, capthick=1.2, markersize=5.5, zorder=10,
        label="Kuchenbaecker et al. (2017) 95% CI"
    )
    ax_c.set_xlabel("Age (years)", fontsize=9.5)
    ax_c.set_ylabel("Cumulative breast cancer risk (%)", fontsize=9.5)
    ax_c.set_xlim(25, 85)
    ax_c.set_ylim(-2, 105)
    ax_c.grid(True, linestyle="--", alpha=0.4)
    ax_c.legend(loc="upper left", fontsize=8.0, frameon=True, facecolor="white", edgecolor="#D1D5DB")

    # Panel d: Normal Tissue Neutral Dynamics
    ax_d = fig.add_subplot(gs[1, 1])
    ax_d.text(-0.16, 1.12, "d", transform=ax_d.transAxes, fontsize=15, fontweight="bold", va="top")
    ax_d.set_title("Normal tissue neutral dynamics", fontsize=11, fontweight="bold", pad=12)

    years = np.linspace(0, 80, 81)
    old_lambda_B = 0.9524
    old_tissue = 1000.0 * (old_lambda_B ** years)

    b_mean = df_stage.groupby("year")["baseline"].mean().values
    norm_mean = (df_stage.groupby("year")["baseline"].mean() +
                 df_stage.groupby("year")["unstable"].mean() +
                 df_stage.groupby("year")["unprotected"].mean()).values

    b_curve = np.insert(b_mean, 0, 1000.0)
    norm_curve = np.insert(norm_mean, 0, 1000.0)

    ax_d.plot(years, norm_curve, color="#059669", lw=2.2, label=r"Calibrated normal cells $N_n(t)$ ($F=1.0$, neutral)")
    ax_d.plot(years, b_curve, color="#2563EB", linestyle="-.", lw=1.8, label=r"Healthy baseline $B(t)$ (retention 85%)")
    ax_d.plot(years, old_tissue, color="#DC2626", linestyle="--", lw=1.8, label=r"Previous model ($d=0.001$, decay $t_{1/2}=14$y)")
    ax_d.axhline(1000, color="#6B7280", linestyle=":", alpha=0.6, label="Initial capacity $N_0 = 1000$")
    ax_d.set_xlabel("Age (years)", fontsize=9.5)
    ax_d.set_ylabel("Normal living cells $N_n(t)$", fontsize=9.5)
    ax_d.set_xlim(0, 80)
    ax_d.set_ylim(0, 1120)
    ax_d.grid(True, linestyle="--", alpha=0.4)
    ax_d.legend(loc="lower left", fontsize=7.8, frameon=True, facecolor="white", edgecolor="#D1D5DB")

    out_pdf = RESULTS_DIR / "Figure3.pdf"
    out_png = RESULTS_DIR / "Figure3.png"
    plt.savefig(out_pdf, bbox_inches="tight")
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Figure 3 guardada en:\n  -> {out_png}")


def render_stage_dynamics():
    df = pd.read_csv(STAGE_CSV)
    grouped = df.groupby("year")
    years = np.array(list(grouped.groups.keys()))

    stages = ["baseline", "unstable", "unprotected", "primer", "tumoral"]
    stage_labels = [
        "Baseline (TP53+/+, BRCA1+/-)",
        "Unstable (TP53+/-, BRCA1+/-)",
        "Unprotected (TP53-/-)",
        "Primer (D1 > 2.0)",
        "Tumoral (Transformed)",
    ]
    stage_colors = ["#2563EB", "#F59E0B", "#10B981", "#8B5CF6", "#DC2626"]
    mean_stages = {s: grouped[s].mean().to_numpy() for s in stages}
    mean_total = grouped["total_alive"].mean().to_numpy()

    fig = plt.figure(figsize=(14.0, 9.5), dpi=300)
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 0.9], hspace=0.32, wspace=0.25)

    # Panel a: Absolute counts
    ax_a = fig.add_subplot(gs[0, 0])
    fig.text(0.04, 0.96, "a", fontsize=15, fontweight="bold", va="top")
    ax_a.set_title("Evolution of Cellular Life Stages over Time (Mean Cell Count, N0=1000)", fontsize=11, fontweight="bold", pad=8)
    y_stack = [mean_stages[s] for s in stages]
    ax_a.stackplot(years, y_stack, labels=stage_labels, colors=stage_colors, alpha=0.88, edgecolor="#111827", linewidth=0.5)
    ax_a.plot(years, mean_total, color="#111827", linestyle="--", linewidth=1.5, label="Total Live Cells")
    ax_a.set_xlim(1, 80)
    ax_a.set_ylim(0, max(mean_total) * 1.15)
    ax_a.set_xlabel("Tissue Age (Years)", fontsize=10)
    ax_a.set_ylabel("Live Cells per Tissue", fontsize=10)
    ax_a.tick_params(labelsize=9)
    ax_a.legend(loc="upper right", fontsize=8.0, frameon=True, edgecolor="#E5E7EB")

    # Panel b: Relative composition
    ax_b = fig.add_subplot(gs[0, 1])
    fig.text(0.52, 0.96, "b", fontsize=15, fontweight="bold", va="top")
    ax_b.set_title("Relative Composition of Cellular Stages (%)", fontsize=11, fontweight="bold", pad=8)
    prop_stack = [mean_stages[s] / np.maximum(mean_total, 1e-9) * 100.0 for s in stages]
    ax_b.stackplot(years, prop_stack, labels=stage_labels, colors=stage_colors, alpha=0.88, edgecolor="#111827", linewidth=0.5)
    ax_b.set_xlim(1, 80)
    ax_b.set_ylim(0, 100)
    ax_b.set_xlabel("Tissue Age (Years)", fontsize=10)
    ax_b.set_ylabel("Percentage of Live Population (%)", fontsize=10)
    ax_b.tick_params(labelsize=9)
    ax_b.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), fontsize=8.0, frameon=False)

    # Panel c: Allelic Dynamics
    ax_c = fig.add_subplot(gs[1, 0])
    fig.text(0.04, 0.48, "c", fontsize=15, fontweight="bold", va="top")
    ax_c.set_title("Genetic Allele Dynamics (BRCA1 & TP53 Breakdown)", fontsize=11, fontweight="bold", pad=8)
    mean_b_pm = grouped["brca1_pm"].mean().to_numpy() / np.maximum(mean_total, 1e-9) * 100.0
    mean_b_mm = grouped["brca1_mm"].mean().to_numpy() / np.maximum(mean_total, 1e-9) * 100.0
    mean_t_pp = grouped["tp53_pp"].mean().to_numpy() / np.maximum(mean_total, 1e-9) * 100.0
    mean_t_pm = grouped["tp53_pm"].mean().to_numpy() / np.maximum(mean_total, 1e-9) * 100.0
    mean_t_mm = grouped["tp53_mm"].mean().to_numpy() / np.maximum(mean_total, 1e-9) * 100.0

    ax_c.plot(years, mean_b_pm, label="BRCA1 +/- (Heterozygous)", color="#2563EB", lw=1.8)
    ax_c.plot(years, mean_b_mm, label="BRCA1 -/- (Null)", color="#DC2626", lw=1.8, linestyle=":")
    ax_c.plot(years, mean_t_pp, label="TP53 +/+ (Wild type)", color="#059669", lw=1.8)
    ax_c.plot(years, mean_t_pm, label="TP53 +/- (First hit)", color="#D97706", lw=1.8, linestyle="--")
    ax_c.plot(years, mean_t_mm, label="TP53 -/- (Double hit)", color="#7C3AED", lw=1.8, linestyle="-.")
    ax_c.set_xlim(1, 80)
    ax_c.set_ylim(0, 105)
    ax_c.set_xlabel("Tissue Age (Years)", fontsize=10)
    ax_c.set_ylabel("Proportion of Population (%)", fontsize=10)
    ax_c.tick_params(labelsize=9)
    ax_c.legend(loc="center right", fontsize=8.0, frameon=True, edgecolor="#E5E7EB")

    # Panel d: Table
    ax_d = fig.add_subplot(gs[1, 1])
    ax_d.axis("off")
    fig.text(0.52, 0.48, "d", fontsize=15, fontweight="bold", va="top")
    ax_d.set_title("Milestone Cell Counts (Mean Cells per Tissue, N0=1000)", fontsize=11, fontweight="bold", pad=8)
    milestones = [20, 40, 60, 80]
    table_data = []
    for m in milestones:
        row = [f"Age {m}"]
        for s in stages:
            row.append(f"{mean_stages[s][m - 1]:.1f}")
        row.append(f"{mean_total[m - 1]:.1f}")
        table_data.append(row)

    col_labels = ["Milestone", "Baseline", "Unstable", "Unprotected", "Primer", "Tumoral", "Total"]
    table = ax_d.table(cellText=table_data, colLabels=col_labels, loc="center", cellLoc="center")
    table.auto_set_font_size(False)
    table.set_fontsize(8.5)
    table.scale(1.05, 1.8)
    for (r, c), cell in table.get_celld().items():
        if r == 0:
            cell.set_facecolor("#E5E7EB")
            cell.set_text_props(fontweight="bold")

    out_png = RESULTS_DIR / "Figure_stage_dynamics.png"
    out_pdf = RESULTS_DIR / "Figure_stage_dynamics.pdf"
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.savefig(out_pdf, bbox_inches="tight")
    plt.close()
    print(f"Figure stage dynamics guardada en:\n  -> {out_png}")


def main():
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    print("Renderizando Figure 4...")
    render_figure4()
    print("Renderizando Figure 3...")
    render_figure3()
    print("Renderizando Figure_stage_dynamics...")
    render_stage_dynamics()
    print("\n[DONE] Todas las figuras a 1000 células generadas con éxito.")


if __name__ == "__main__":
    main()
