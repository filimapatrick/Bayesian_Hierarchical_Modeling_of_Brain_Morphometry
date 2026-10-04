#!/usr/bin/env python3
"""
Experiment 2: Real-World Clinical Measurement Validity & Failure Boundaries

Maps the operational feasibility boundaries of automated macro-morphometry across
the heterogeneous Nigerian clinical dataset (N=218 scans).

Investigates:
1. Pipeline retention rate across clinical cohorts vs. standard micro-segmentation (FreeSurfer pilot benchmark).
2. Logistic failure probability: P(Pipeline Unviability / Failure) = f(SliceThickness, Contrast, FieldStrength).
3. Evidence-based feasibility boundaries for opportunistic LMIC clinical archives.
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path
from typing import Dict, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats

PROJECT_ROOT = Path(__file__).resolve().parent.parent

# Ensure scipy.signal.gaussian compatibility for ArviZ under SciPy 1.13+
import scipy.signal
import scipy.signal.windows
if not hasattr(scipy.signal, "gaussian"):
    scipy.signal.gaussian = scipy.signal.windows.gaussian

import arviz as az
import pymc as pm

def fit_logistic_failure_model(
    df: pd.DataFrame,
    draws: int = 1500,
    tune: int = 1000,
    chains: int = 4,
    random_seed: int = 42,
) -> Tuple[pd.DataFrame, Dict[str, float], Dict[str, np.ndarray]]:
    """
    Fits a Bayesian multivariate logistic regression model predicting pipeline unviability / measurement failure
    as a function of slice thickness, gadolinium contrast, and low field strength using PyMC NUTS.

    Likelihood:
      y_i ~ Bernoulli(p_i)
      logit(p_i) = alpha + beta_thick * (h_i - 1.0) + beta_contrast * Contrast_i + beta_lowfield * LowField_i

    Prespecified Quality Control Failure Criteria:
      1. Scout / localizer series (< 10 slices, truncated FOV)
      2. Evans' Index: < 0.18 or > 0.85 (Anatomically implausible automated geometry)
      3. PEF: < 0.45 or > 0.96 (Envelope extraction boundary violations)
      4. VBR: <= 0.0 or > 0.60 (Non-physical / extreme segmentation leak)
      5. Hemispheric Asymmetry Index (HAI): > 35.0% (Gross unilateral coil cutoff or tilt)
    """
    # 1. Scout / localizer series
    scout_fail = (df["is_volumetric_valid"] == 0.0)
    
    # 2. Biological boundary violations & extraction failures
    pef_col = "pef" if "pef" in df.columns else "bpf"
    ei_fail = df["evans_index"].isna() | (df["evans_index"] < 0.18) | (df["evans_index"] > 0.85)
    pef_fail = df[pef_col].isna() | (df[pef_col] < 0.45) | (df[pef_col] > 0.96)
    vbr_fail = df["vbr"].isna() | (df["vbr"] <= 0.0) | (df["vbr"] > 0.60)
    asym_fail = df["asymmetry_index"].isna() | (df["asymmetry_index"] > 35.0)
    
    # Combined measurement unviability indicator
    is_failed = scout_fail | ei_fail | pef_fail | vbr_fail | asym_fail
    y = is_failed.astype(int).values
    
    # Covariates (centered thickness at 1.0 mm reference)
    x_thick = df["slice_thickness_mm"].fillna(df["slice_thickness_mm"].median()).values
    x_contrast = df["contrast_enhanced"].astype(float).values
    x_lowfield = (df["magnetic_field_strength"] <= 0.35).astype(float).values
    
    with pm.Model() as model:
        alpha = pm.Normal("alpha", mu=-2.0, sigma=2.0)
        beta_thick = pm.Normal("beta_thick", mu=0.0, sigma=1.0)
        beta_contrast = pm.Normal("beta_contrast", mu=0.0, sigma=1.0)
        beta_lowfield = pm.Normal("beta_lowfield", mu=0.0, sigma=1.0)
        
        logit_p = alpha + beta_thick * (x_thick - 1.0) + beta_contrast * x_contrast + beta_lowfield * x_lowfield
        pm.Bernoulli("obs", logit_p=logit_p, observed=y)
        
        idata = pm.sample(
            draws=draws,
            tune=tune,
            chains=chains,
            target_accept=0.95,
            random_seed=random_seed,
            progressbar=False,
        )
        
    summary = az.summary(idata, hdi_prob=0.95)
    
    post_alpha = idata.posterior["alpha"].values.flatten()
    post_thick = idata.posterior["beta_thick"].values.flatten()
    post_contrast = idata.posterior["beta_contrast"].values.flatten()
    post_lowfield = idata.posterior["beta_lowfield"].values.flatten()
    
    or_thick = np.exp(post_thick)
    or_contrast = np.exp(post_contrast)
    or_lowfield = np.exp(post_lowfield)
    
    records = [
        {
            "Predictor": "Intercept (1.0mm, Unenhanced, 1.5T)",
            "Mean_LogOdds": float(summary.loc["alpha", "mean"]),
            "SD": float(summary.loc["alpha", "sd"]),
            "HDI_2.5%": float(summary.loc["alpha", "hdi_2.5%"]),
            "HDI_97.5%": float(summary.loc["alpha", "hdi_97.5%"]),
            "Odds_Ratio_Median": float(np.median(np.exp(post_alpha))),
            "OR_HDI_2.5%": float(np.percentile(np.exp(post_alpha), 2.5)),
            "OR_HDI_97.5%": float(np.percentile(np.exp(post_alpha), 97.5)),
            "R_hat": float(summary.loc["alpha", "r_hat"]),
            "ESS_bulk": float(summary.loc["alpha", "ess_bulk"]),
        },
        {
            "Predictor": "Slice Thickness (per mm above 1mm)",
            "Mean_LogOdds": float(summary.loc["beta_thick", "mean"]),
            "SD": float(summary.loc["beta_thick", "sd"]),
            "HDI_2.5%": float(summary.loc["beta_thick", "hdi_2.5%"]),
            "HDI_97.5%": float(summary.loc["beta_thick", "hdi_97.5%"]),
            "Odds_Ratio_Median": float(np.median(or_thick)),
            "OR_HDI_2.5%": float(np.percentile(or_thick, 2.5)),
            "OR_HDI_97.5%": float(np.percentile(or_thick, 97.5)),
            "R_hat": float(summary.loc["beta_thick", "r_hat"]),
            "ESS_bulk": float(summary.loc["beta_thick", "ess_bulk"]),
        },
        {
            "Predictor": "Contrast (+C)",
            "Mean_LogOdds": float(summary.loc["beta_contrast", "mean"]),
            "SD": float(summary.loc["beta_contrast", "sd"]),
            "HDI_2.5%": float(summary.loc["beta_contrast", "hdi_2.5%"]),
            "HDI_97.5%": float(summary.loc["beta_contrast", "hdi_97.5%"]),
            "Odds_Ratio_Median": float(np.median(or_contrast)),
            "OR_HDI_2.5%": float(np.percentile(or_contrast, 2.5)),
            "OR_HDI_97.5%": float(np.percentile(or_contrast, 97.5)),
            "R_hat": float(summary.loc["beta_contrast", "r_hat"]),
            "ESS_bulk": float(summary.loc["beta_contrast", "ess_bulk"]),
        },
        {
            "Predictor": "Low Field (<=0.35T)",
            "Mean_LogOdds": float(summary.loc["beta_lowfield", "mean"]),
            "SD": float(summary.loc["beta_lowfield", "sd"]),
            "HDI_2.5%": float(summary.loc["beta_lowfield", "hdi_2.5%"]),
            "HDI_97.5%": float(summary.loc["beta_lowfield", "hdi_97.5%"]),
            "Odds_Ratio_Median": float(np.median(or_lowfield)),
            "OR_HDI_2.5%": float(np.percentile(or_lowfield, 2.5)),
            "OR_HDI_97.5%": float(np.percentile(or_lowfield, 97.5)),
            "R_hat": float(summary.loc["beta_lowfield", "r_hat"]),
            "ESS_bulk": float(summary.loc["beta_lowfield", "ess_bulk"]),
        },
    ]
    summary_df = pd.DataFrame(records)
    
    model_params = {
        "alpha": float(summary.loc["alpha", "mean"]),
        "beta_thick": float(summary.loc["beta_thick", "mean"]),
        "beta_contrast": float(summary.loc["beta_contrast", "mean"]),
        "beta_lowfield": float(summary.loc["beta_lowfield", "mean"]),
    }
    traces = {
        "alpha": post_alpha,
        "beta_thick": post_thick,
        "beta_contrast": post_contrast,
        "beta_lowfield": post_lowfield,
    }
    return summary_df, model_params, traces


def generate_feasibility_plot(
    df: pd.DataFrame,
    logit_summary: pd.DataFrame,
    model_params: Dict[str, float],
    traces: Dict[str, np.ndarray],
    plot_path: Path,
):
    """
    Renders publication Figure 2:
    Clinical Feasibility Boundaries, Failure Rates, and Pipeline Retention Contrast.
    """
    plot_path.parent.mkdir(parents=True, exist_ok=True)
    sns.set_theme(style="whitegrid", font_scale=1.05)
    
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    
    # -------------------------------------------------------------
    # Panel A: Pipeline Retention by Diagnostic Cohort (Macro vs FreeSurfer Pilot)
    # -------------------------------------------------------------
    ax1 = axes[0, 0]
    cohorts = ["Control", "Dementia", "Epilepsy", "Hydrocephalus", "Parkinson"]
    
    # Macro-Morphometry Engine retention
    macro_rates = []
    for c in cohorts:
        sub_df = df[df["diagnosis"] == c.upper()]
        rate = (sub_df["is_volumetric_valid"] == 1.0).mean() * 100.0 if len(sub_df) > 0 else 0.0
        macro_rates.append(rate)
        
    # FreeSurfer Pilot Attrition Benchmark (calculated from results/tables/freesurfer_qc.csv)
    fs_qc_path = PROJECT_ROOT / "results" / "tables" / "freesurfer_qc.csv"
    if not fs_qc_path.exists():
        raise FileNotFoundError(
            f"FreeSurfer QC benchmark file not found at {fs_qc_path}. "
            "Scientific reproducibility requires participant-level QC data; hard-coded fallbacks are disabled."
        )
    fs_df = pd.read_csv(fs_qc_path)
    fs_rates = []
    for c in cohorts:
        sub_fs = fs_df[fs_df["diagnosis"] == c.upper()]
        rate = (sub_fs["qc_pass"] == True).mean() * 100.0 if len(sub_fs) > 0 else 0.0
        fs_rates.append(rate)
    
    x = np.arange(len(cohorts))
    width = 0.36
    
    b1 = ax1.bar(x - width/2, macro_rates, width, label="Macro-Morphometry Engine", color="#1f77b4", alpha=0.9, edgecolor="black", linewidth=1.2)
    b2 = ax1.bar(x + width/2, fs_rates, width, label="FreeSurfer (Pilot Benchmark)", color="#d62728", alpha=0.8, edgecolor="black", linewidth=1.2)
    
    ax1.set_ylabel("Cohort Retention Rate (%)", fontweight="bold")
    ax1.set_title("A. Pipeline Retention Rate Across Clinical Cohorts", fontweight="bold", pad=12)
    ax1.set_xticks(x)
    ax1.set_xticklabels(cohorts, fontweight="bold")
    ax1.set_ylim(0, 115)
    ax1.axhline(100, color="gray", linestyle="--", alpha=0.5)
    ax1.legend(loc="upper right", frameon=True, framealpha=0.9)
    
    for bar in b1:
        h = bar.get_height()
        ax1.annotate(f"{h:.1f}%", xy=(bar.get_x() + bar.get_width() / 2, h),
                     xytext=(0, 4), textcoords="offset points", ha="center", va="bottom", fontsize=8.5, fontweight="bold")
    for bar in b2:
        h = bar.get_height()
        ax1.annotate(f"{h:.1f}%", xy=(bar.get_x() + bar.get_width() / 2, h),
                     xytext=(0, 4), textcoords="offset points", ha="center", va="bottom", fontsize=8.5, fontweight="bold", color="#a00")
        
    # -------------------------------------------------------------
    # Panel B: Modeled Probability of Measurement Failure vs Slice Thickness
    # -------------------------------------------------------------
    ax2 = axes[0, 1]
    th_range = np.linspace(1.0, 8.0, 150)
    delta_th = th_range[None, :] - 1.0  # shape: (1, 150)
    
    alpha_samples = traces["alpha"][:, None]          # shape: (S, 1)
    b_thick_samples = traces["beta_thick"][:, None]   # shape: (S, 1)
    b_cont_samples = traces["beta_contrast"][:, None] # shape: (S, 1)
    b_low_samples = traces["beta_lowfield"][:, None]  # shape: (S, 1)
    
    # 1. Unenhanced, 1.5T
    logit_p1_samples = alpha_samples + b_thick_samples * delta_th
    p1_samples = 1.0 / (1.0 + np.exp(-logit_p1_samples))
    p1_median = np.median(p1_samples, axis=0)
    p1_lo = np.percentile(p1_samples, 2.5, axis=0)
    p1_hi = np.percentile(p1_samples, 97.5, axis=0)
    
    # 2. Contrast-enhanced (+C), 1.5T
    logit_p2_samples = alpha_samples + b_thick_samples * delta_th + b_cont_samples
    p2_samples = 1.0 / (1.0 + np.exp(-logit_p2_samples))
    p2_median = np.median(p2_samples, axis=0)
    p2_lo = np.percentile(p2_samples, 2.5, axis=0)
    p2_hi = np.percentile(p2_samples, 97.5, axis=0)
    
    # 3. Low-field (0.35T), Unenhanced
    logit_p3_samples = alpha_samples + b_thick_samples * delta_th + b_low_samples
    p3_samples = 1.0 / (1.0 + np.exp(-logit_p3_samples))
    p3_median = np.median(p3_samples, axis=0)
    p3_lo = np.percentile(p3_samples, 2.5, axis=0)
    p3_hi = np.percentile(p3_samples, 97.5, axis=0)
    
    ax2.plot(th_range, p1_median, label="1.5T Unenhanced (Median)", color="#2ca02c", linewidth=2.5)
    ax2.fill_between(th_range, p1_lo, p1_hi, color="#2ca02c", alpha=0.18, label="1.5T Unenhanced 95% HDI")
    
    ax2.plot(th_range, p2_median, label="1.5T Contrast (+C) (Median)", color="#ff7f0e", linewidth=2.5, linestyle="--")
    ax2.fill_between(th_range, p2_lo, p2_hi, color="#ff7f0e", alpha=0.18, label="1.5T Contrast 95% HDI")
    
    ax2.plot(th_range, p3_median, label="0.35T Low-Field (Median)", color="#9467bd", linewidth=2.5, linestyle="-.")
    ax2.fill_between(th_range, p3_lo, p3_hi, color="#9467bd", alpha=0.18, label="0.35T Low-Field 95% HDI")
    
    ax2.set_xlabel("Acquired Slice Thickness (mm)", fontweight="bold")
    ax2.set_ylabel("P(Measurement Failure / Unviable)", fontweight="bold")
    ax2.set_title("B. Modeled Failure Probability with 95% Posterior HDI Bands", fontweight="bold", pad=12)
    ax2.set_xlim(1.0, 8.0)
    ax2.set_ylim(-0.02, 1.02)
    ax2.axvline(5.0, color="red", linestyle=":", alpha=0.7, label="Standard Clinical 2D (5mm)")
    ax2.legend(loc="upper left", frameon=True, framealpha=0.9, fontsize=8.5)
    
    # -------------------------------------------------------------
    # Panel C: Log-Odds of Measurement Failure Risk Factors
    # -------------------------------------------------------------
    ax3 = axes[1, 0]
    preds = logit_summary.iloc[1:].copy()
    y_pos = np.arange(len(preds))
    
    means = preds["Mean_LogOdds"].values
    hdi_lo = preds["HDI_2.5%"].values
    hdi_hi = preds["HDI_97.5%"].values
    xerr = [means - hdi_lo, hdi_hi - means]
    
    ax3.errorbar(means, y_pos, xerr=xerr, fmt='o', color="#1f77b4",
                 ecolor="#333333", elinewidth=2, capsize=5, markersize=8)
    
    ax3.axvline(0.0, color="red", linestyle="--", alpha=0.7)
    ax3.set_yticks(y_pos)
    ax3.set_yticklabels(preds["Predictor"].values, fontweight="bold")
    ax3.set_xlabel("Posterior Log-Odds [95% HDI]", fontweight="bold")
    ax3.set_title("C. Bayesian Adjusted Log-Odds for Structural Failure", fontweight="bold", pad=12)
    
    for idx, (_, row) in enumerate(preds.iterrows()):
        ax3.annotate(f"OR = {row['Odds_Ratio_Median']:.2f}\n[95% HDI: {row['OR_HDI_2.5%']:.2f}, {row['OR_HDI_97.5%']:.2f}]",
                     xy=(row['Mean_LogOdds'], idx), xytext=(12, -8),
                     textcoords="offset points", fontsize=8.0, fontweight="bold")

    # -------------------------------------------------------------
    # Panel D: Feasibility Decision Matrix Table
    # -------------------------------------------------------------
    ax4 = axes[1, 1]
    ax4.axis("off")
    
    # Posterior predictions for representative clinical acquisition profiles:
    # 1. 3D High-Res (1.0 mm, unenhanced, 1.5T)
    prof1_p = 1.0 / (1.0 + np.exp(-(traces["alpha"])))
    prof1_med = np.median(prof1_p) * 100.0
    prof1_lo = np.percentile(prof1_p, 2.5) * 100.0
    prof1_hi = np.percentile(prof1_p, 97.5) * 100.0

    # 2. 2D Standard (4.0 mm, unenhanced, 1.5T)
    prof2_p = 1.0 / (1.0 + np.exp(-(traces["alpha"] + traces["beta_thick"] * 3.0)))
    prof2_med = np.median(prof2_p) * 100.0
    prof2_lo = np.percentile(prof2_p, 2.5) * 100.0
    prof2_hi = np.percentile(prof2_p, 97.5) * 100.0

    # 3. 2D Thick-Slice (5.0 mm, +C, 1.5T)
    prof3_p = 1.0 / (1.0 + np.exp(-(traces["alpha"] + traces["beta_thick"] * 4.0 + traces["beta_contrast"])))
    prof3_med = np.median(prof3_p) * 100.0
    prof3_lo = np.percentile(prof3_p, 2.5) * 100.0
    prof3_hi = np.percentile(prof3_p, 97.5) * 100.0

    # 4. Low-Field (5.0 mm, unenhanced, 0.35T)
    prof4_p = 1.0 / (1.0 + np.exp(-(traces["alpha"] + traces["beta_thick"] * 4.0 + traces["beta_lowfield"])))
    prof4_med = np.median(prof4_p) * 100.0
    prof4_lo = np.percentile(prof4_p, 2.5) * 100.0
    prof4_hi = np.percentile(prof4_p, 97.5) * 100.0

    # FreeSurfer empirical benchmark
    fs_fail_rate = (1.0 - (fs_df["qc_pass"] == True).mean()) * 100.0
    
    matrix_data = [
        ["Acquisition Profile", "Feasible Biomarkers", "Posterior Failure Risk\n[Median, 95% HDI]", "Recommended\nStatistical Handling"],
        ["3D High-Res\n(1.0 mm, 1.5T)", "Evans' Index, PEF,\nVBR, Subcortical", f"{prof1_med:.2f}%\n[{prof1_lo:.2f}%, {prof1_hi:.2f}%]", "Direct pooling"],
        ["2D Standard\n(4.0 mm, 1.5T)", "Evans' Index (Flagship),\nPEF", f"{prof2_med:.2f}%\n[{prof2_lo:.2f}%, {prof2_hi:.2f}%]", "Heteroskedastic noise\nscaling (λ)"],
        ["2D Thick-Slice\n(5.0 mm, +C, 1.5T)", "Evans' Index only", f"{prof3_med:.2f}%\n[{prof3_lo:.2f}%, {prof3_hi:.2f}%]", "Contrast covariate (δ)\n+ noise scaling (λ)"],
        ["Low-Field\n(5.0 mm, 0.35T)", "Evans' Index\n(with manual QC)", f"{prof4_med:.2f}%\n[{prof4_lo:.2f}%, {prof4_hi:.2f}%]", "Institutional random\neffect + QC bounds"],
        ["FreeSurfer Pilot\n(2D Thick-Slice)", "None\n(Catastrophic attrition)", f"{fs_fail_rate:.1f}%\n(Empirical failure)", "DO NOT DEPLOY\n(Biologically unviable)"]
    ]
    
    tbl = ax4.table(cellText=matrix_data, loc="center", cellLoc="center",
                    colWidths=[0.24, 0.26, 0.25, 0.25])
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(8.2)
    tbl.scale(1.0, 2.2)
    
    for (row, col), cell in tbl.get_celld().items():
        cell.set_edgecolor("#cccccc")
        if row == 0:
            cell.set_facecolor("#2c3e50")
            cell.get_text().set_color("white")
            cell.get_text().set_weight("bold")
        else:
            if row == 5:
                cell.set_facecolor("#ffebee")
                cell.get_text().set_color("#c62828")
                cell.get_text().set_weight("bold")
            else:
                cell.set_facecolor("#fafafa" if row % 2 == 0 else "#ffffff")
                
    ax4.set_title("D. Opportunistic Clinical Neuroimaging Feasibility Matrix", fontweight="bold", pad=12)

    plt.savefig(plot_path, dpi=300, bbox_inches="tight")
    plt.close()


def main():
    parser = argparse.ArgumentParser(description="Run Experiment 2: Clinical Feasibility & Failure Boundaries.")
    parser.add_argument("--features_csv", type=str, default="results/tables/macro_features.csv")
    parser.add_argument("--output_csv", type=str, default="results/tables/experiment2_failure_boundaries.csv")
    parser.add_argument("--plot_path", type=str, default="results/figures/figure2_clinical_feasibility_boundaries.png")
    args = parser.parse_args()

    feat_path = Path(args.features_csv)
    out_csv = Path(args.output_csv)
    plot_path = Path(args.plot_path)

    if not feat_path.exists():
        sys.exit(f"Error: {feat_path} does not exist. Run feature extraction first.")

    t_start = time.time()
    df = pd.read_csv(feat_path)

    print("=" * 75)
    print("EXPERIMENT 2: REAL-WORLD CLINICAL FEASIBILITY & FAILURE BOUNDARIES")
    print(f"Dataset: {len(df)} clinical scans from {feat_path}")
    print("=" * 75)

    total_scans = len(df)
    valid_scans = (df["is_volumetric_valid"] == 1.0).sum()
    failed_scans = (df["is_volumetric_valid"] == 0.0).sum()
    retention_rate = (valid_scans / total_scans) * 100.0

    print(f"\n📊 Global Pipeline Retention:")
    print(f"  • Total Archive Records: {total_scans}")
    print(f"  • Valid Analyzable Scans: {valid_scans} ({retention_rate:.1f}%)")
    print(f"  • Flagged Scouts / Invalids: {failed_scans} ({100.0 - retention_rate:.1f}%)")

    logit_summary, model_params, traces = fit_logistic_failure_model(df)
    print("\n📈 Bayesian Logistic Regression Failure Model [Logit(P(Failure))]:")
    print(logit_summary.to_string(index=False))

    out_csv.parent.mkdir(parents=True, exist_ok=True)
    logit_summary.to_csv(out_csv, index=False)
    print(f"\n✓ Saved failure model table: {out_csv}")

    generate_feasibility_plot(df, logit_summary, model_params, traces, plot_path)
    print(f"✓ Saved publication figure: {plot_path}")
    print(f"Experiment 2 completed in {time.time() - t_start:.2f}s")
    print("=" * 75)


if __name__ == "__main__":
    main()
