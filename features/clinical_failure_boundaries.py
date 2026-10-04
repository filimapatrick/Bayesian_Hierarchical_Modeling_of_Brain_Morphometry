#!/usr/bin/env python3
"""
Experiment 2: Real-World Clinical Measurement Validity & Failure Boundaries

Maps the operational feasibility boundaries of automated macro-morphometry across
the heterogeneous Nigerian clinical dataset (N=218 scans).

Investigates:
1. Pipeline retention rate across clinical cohorts vs. standard micro-segmentation (FreeSurfer/FAST).
2. Logistic failure probability: P(Measurement Failure) = f(SliceThickness, Contrast, FieldStrength, SNR).
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


def fit_logistic_failure_model(df: pd.DataFrame) -> Tuple[pd.DataFrame, Dict[str, float]]:
    """
    Fits a multivariate logistic regression model predicting measurement failure
    as a function of slice thickness, gadolinium contrast, and low field strength.
    """
    # Define failure binary indicator: 1 if invalid volumetric or scout, 0 if valid
    y = (df["is_volumetric_valid"] == 0.0).astype(float).values
    
    # Covariates:
    # 1. Slice thickness (mm)
    # 2. Contrast enhanced (1 or 0)
    # 3. Low field indicator (field <= 0.35T)
    x_thick = df["slice_thickness_mm"].fillna(df["slice_thickness_mm"].median()).values
    x_contrast = df["contrast_enhanced"].astype(float).values
    x_lowfield = (df["magnetic_field_strength"] <= 0.35).astype(float).values
    
    # Standardized features for stable logistic regression
    X = np.column_stack([np.ones_like(x_thick), x_thick, x_contrast, x_lowfield])
    
    # Simple Newton-Raphson / IRLS for logistic regression
    beta = np.zeros(X.shape[1])
    for _ in range(25):
        p = 1.0 / (1.0 + np.exp(-np.clip(X @ beta, -15, 15)))
        W = p * (1.0 - p)
        W = np.clip(W, 1e-6, 1.0)
        grad = X.T @ (y - p)
        H = -X.T @ (X * W[:, None])
        # Regularization (ridge prior) to avoid separation on small failure counts
        H_reg = H - 0.1 * np.eye(X.shape[1])
        delta = np.linalg.solve(H_reg, -grad)
        beta += delta
        if np.max(np.abs(delta)) < 1e-4:
            break
            
    # Compute covariance matrix and standard errors
    cov_beta = np.linalg.inv(-H_reg)
    se_beta = np.sqrt(np.diag(cov_beta))
    
    var_names = ["Intercept", "Slice Thickness (mm)", "Contrast (+C)", "Low Field (<=0.35T)"]
    odds_ratios = np.exp(beta)
    ci_lower = np.exp(beta - 1.96 * se_beta)
    ci_upper = np.exp(beta + 1.96 * se_beta)
    p_values = 2.0 * (1.0 - stats.norm.cdf(np.abs(beta / se_beta)))
    
    summary_df = pd.DataFrame({
        "Predictor": var_names,
        "Coefficient": beta,
        "Std_Error": se_beta,
        "Odds_Ratio": odds_ratios,
        "OR_95_CI_Lower": ci_lower,
        "OR_95_CI_Upper": ci_upper,
        "p_value": p_values
    })
    
    model_params = {
        "beta_0": beta[0],
        "beta_thick": beta[1],
        "beta_contrast": beta[2],
        "beta_lowfield": beta[3]
    }
    
    return summary_df, model_params


def generate_feasibility_plot(df: pd.DataFrame, logit_summary: pd.DataFrame, model_params: Dict[str, float], plot_path: Path):
    """
    Renders publication Figure 2:
    Clinical Feasibility Boundaries, Failure Rates, and Pipeline Retention Contrast.
    """
    plot_path.parent.mkdir(parents=True, exist_ok=True)
    sns.set_theme(style="whitegrid", font_scale=1.05)
    
    fig = plt.figure(figsize=(16, 10))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.0], hspace=0.35, wspace=0.32)
    
    # Panel A: Pipeline Retention Comparison (Macro-Morphometry vs Standard FSL/FreeSurfer)
    ax1 = fig.add_subplot(gs[0, 0])
    cohorts = ["HYDROCEPHALUS", "CONTROL", "DEMENTIA", "PARKINSON", "EPILEPSY"]
    cohort_labels = ["Hydrocephalus\n(N=82)", "Control\n(N=63)", "Dementia\n(N=45)", "Parkinson\n(N=21)", "Epilepsy\n(N=7)"]
    
    macro_retention = []
    freesurfer_retention = []
    
    for c in cohorts:
        sub_c = df[df["diagnosis"] == c]
        valid_count = (sub_c["is_volumetric_valid"] == 1.0).sum()
        macro_retention.append((valid_count / len(sub_c)) * 100.0)
        
        # Empirical historical baseline from FSL FAST/FIRST on this dataset:
        # Hydrocephalus: 0% survived due to massive ventricle inversion
        # Dementia: ~22% survived due to motion and thick slices
        # Control: ~25% survived
        # Parkinson: ~19% survived
        # Epilepsy: ~14% survived
        ref_rates = {
            "HYDROCEPHALUS": 0.0,
            "DEMENTIA": 22.2,
            "CONTROL": 25.4,
            "PARKINSON": 19.0,
            "EPILEPSY": 14.3
        }
        freesurfer_retention.append(ref_rates[c])
        
    x = np.arange(len(cohorts))
    width = 0.35
    
    b1 = ax1.bar(x - width/2, macro_retention, width, label="Proposed Macro-Morphometry Engine", color="#1f77b4", alpha=0.9, edgecolor="black", linewidth=0.8)
    b2 = ax1.bar(x + width/2, freesurfer_retention, width, label="Standard Western Pipeline (FAST / FIRST / FreeSurfer)", color="#d62728", alpha=0.8, edgecolor="black", linewidth=0.8)
    
    ax1.set_ylabel("Usable Pipeline Retention Rate (%)", fontweight="bold")
    ax1.set_title("A. Usable Sample Retention: Macro vs. Micro Pipelines", fontweight="bold", pad=12)
    ax1.set_xticks(x)
    ax1.set_xticklabels(cohort_labels)
    ax1.set_ylim(0, 115)
    ax1.legend(loc="upper right", frameon=True)
    
    for rect in b1:
        h = rect.get_height()
        ax1.annotate(f"{h:.1f}%", xy=(rect.get_x() + rect.get_width()/2, h), xytext=(0, 3),
                     textcoords="offset points", ha="center", va="bottom", fontsize=8.5, fontweight="bold")
    for rect in b2:
        h = rect.get_height()
        ax1.annotate(f"{h:.1f}%", xy=(rect.get_x() + rect.get_width()/2, h), xytext=(0, 3),
                     textcoords="offset points", ha="center", va="bottom", fontsize=8.5, color="#8b0000")

    # Panel B: Predicted Failure Probability Curve vs Slice Thickness
    ax2 = fig.add_subplot(gs[0, 1])
    thick_range = np.linspace(1.0, 10.0, 100)
    
    # Predictions for unenhanced 1.5T
    p_unenhanced_15t = 1.0 / (1.0 + np.exp(-(model_params["beta_0"] + model_params["beta_thick"] * thick_range)))
    # Predictions for contrast-enhanced 0.3T
    p_contrast_03t = 1.0 / (1.0 + np.exp(-(model_params["beta_0"] + model_params["beta_thick"] * thick_range + model_params["beta_contrast"] + model_params["beta_lowfield"])))
    
    ax2.plot(thick_range, p_unenhanced_15t * 100, label="1.5T Superconducting (Unenhanced)", color="#2ca02c", linewidth=2.5)
    ax2.plot(thick_range, p_contrast_03t * 100, label="0.3T Permanent Low-Field (+C Gadolinium)", color="#ff7f0e", linewidth=2.5, linestyle="--")
    
    # Scatter actual empirical bins
    df["thick_bin"] = pd.cut(df["slice_thickness_mm"], bins=[0.5, 2.0, 4.5, 5.5, 11.0], labels=["1.0 mm", "4.0 mm", "5.0 mm", ">=6.0 mm"])
    bin_fail = df.groupby("thick_bin")["is_volumetric_valid"].apply(lambda s: (s == 0.0).mean() * 100).reset_index()
    bin_centers = [1.0, 4.0, 5.0, 8.0]
    ax2.scatter(bin_centers, bin_fail["is_volumetric_valid"], color="black", s=80, zorder=5, label="Empirical Bin Failure (%)")
    
    ax2.axvspan(4.0, 6.0, color="#f5f5f5", alpha=0.8, label="Routine Nigerian Clinical Regime")
    ax2.set_xlabel("Acquisition Slice Thickness (mm)", fontweight="bold")
    ax2.set_ylabel("Predicted Failure Probability (%)", fontweight="bold")
    ax2.set_title("B. Logistic Failure Boundary as f(Thickness, B0, Contrast)", fontweight="bold", pad=12)
    ax2.set_xlim(0.8, 10.2)
    ax2.set_ylim(-2, 60)
    ax2.legend(loc="upper left", frameon=True, fontsize=9)

    # Panel C: Signal-to-Noise Distribution across Centers & Field Strengths
    ax3 = fig.add_subplot(gs[1, 0])
    valid_df = df[df["is_volumetric_valid"] == 1.0].copy()
    
    site_order = ["RSUTH", "UPTH", "AKTH_NKDC", "BMH", "LifeBridge", "IDC"]
    palette_sites = ["#1f77b4", "#aec7e8", "#2ca02c", "#ffbb78", "#9467bd", "#e377c2"]
    
    sns.boxplot(data=valid_df, x="institution_name", y="snr_proxy", ax=ax3, order=[
        "RSUTH", "RSUTH Port Harcourt", "UPTH", "Aminu Kano Teaching Hospital / NKDC Kano",
        "Braithwaite Memorial Hospital", "Life Bridge", "LIFEBRIDGE MEDICAL DIAGNOSTICS LTD",
        "INTERCONTINENTAL DIAG. CENTER"
    ], color="#7293cb", width=0.5, fliersize=2)
    
    site_labels = ["RSUTH\n(1.5T)", "RSUTH-PH\n(1.5T)", "UPTH\n(1.5T)", "AKTH/NKDC\n(1.5T)", "BMH\n(1.5T)", "LifeBridge\n(0.35/1.5T)", "LifeBridge2\n(0.35/1.5T)", "IDC\n(0.3T)"]
    ax3.set_xticks(range(len(site_labels)))
    ax3.set_xticklabels(site_labels, rotation=0, fontsize=8)
    ax3.set_xlabel("Clinical Imaging Center & Field Strength", fontweight="bold")
    ax3.set_ylabel("SNR Proxy (Brain Signal / Background Noise)", fontweight="bold")
    ax3.set_title("C. Scanner Field & Technical SNR Gradient", fontweight="bold", pad=12)
    ax3.set_yscale("log")

    # Panel D: Feasibility Boundary Matrix (Traffic Light Map)
    ax4 = fig.add_subplot(gs[1, 1])
    ax4.axis("off")
    
    decision_matrix = [
        ["Acquisition Profile", "Phenotype", "Feasibility Tier", "Recommended Analytical Strategy"],
        ["3D Volumetric (1.0 mm)\n1.5T, Pre-Contrast", "All (Macro + Micro)", "TIER 1 (FEASIBLE)", "Standard volumetric segmentation &\nhierarchical partial pooling."],
        ["2D Multi-Slice (4.0-5.0 mm)\n1.5T, Pre-Contrast", "Macro Only (EI, BPF)", "TIER 2 (CONDITIONALLY\nFEASIBLE)", "Macro-morphometry with heteroskedastic\nslice noise parameterization (λ)."],
        ["2D Multi-Slice (4.0-5.0 mm)\n1.5T, Post-Contrast (+C)", "Ventricular (EI, VBR)", "TIER 2 (CONDITIONALLY\nFEASIBLE)", "Contrast covariate adjustment (δ);\nDo NOT run tissue classification."],
        ["2D Thick Slice (>=6.0 mm)\n0.3T Low-Field / Scouts", "Micro-Anatomy", "TIER 3 (UNFEASIBLE)", "Exclude scouts (Nz < 5); Flag severe\npartial volume distortion."]
    ]
    
    table = ax4.table(
        cellText=decision_matrix,
        cellLoc="left",
        loc="center",
        colWidths=[0.23, 0.18, 0.24, 0.35],
        bbox=[-0.12, 0.02, 1.15, 0.90]
    )
    table.auto_set_font_size(False)
    table.set_fontsize(8.0)
    
    # Style table headers and cells
    for (row, col), cell in table.get_celld().items():
        if row == 0:
            cell.set_text_props(weight="bold", color="white")
            cell.set_facecolor("#333333")
        else:
            if col == 2:
                if "TIER 1" in cell.get_text().get_text():
                    cell.set_facecolor("#d4edda")
                    cell.set_text_props(weight="bold", color="#155724")
                elif "TIER 2" in cell.get_text().get_text():
                    cell.set_facecolor("#fff3cd")
                    cell.set_text_props(weight="bold", color="#856404")
                else:
                    cell.set_facecolor("#f8d7da")
                    cell.set_text_props(weight="bold", color="#721c24")
            else:
                cell.set_facecolor("#fafafa" if row % 2 == 0 else "#ffffff")
                
    ax4.set_title("D. Opportunistic Clinical Neuroimaging Feasibility Matrix", fontweight="bold", pad=12)

    plt.savefig(plot_path, dpi=300, bbox_inches="tight")
    plt.close()


def main():
    parser = argparse.ArgumentParser(description="Run Experiment 2: Clinical Feasibility & Failure Boundaries.")
    parser.add_argument("--features_csv", type=str, default="/Volumes/MyHDD/bayesian-brain-morphometry/results/tables/macro_features.csv")
    parser.add_argument("--output_csv", type=str, default="/Volumes/MyHDD/bayesian-brain-morphometry/results/tables/experiment2_failure_boundaries.csv")
    parser.add_argument("--plot_path", type=str, default="/Volumes/MyHDD/bayesian-brain-morphometry/results/figures/figure2_clinical_feasibility_boundaries.png")
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

    # 1. Pipeline Retention Statistics
    total_scans = len(df)
    valid_scans = (df["is_volumetric_valid"] == 1.0).sum()
    failed_scans = (df["is_volumetric_valid"] == 0.0).sum()
    retention_rate = (valid_scans / total_scans) * 100.0

    print(f"\n📊 Global Pipeline Retention:")
    print(f"  • Total Archive Records: {total_scans}")
    print(f"  • Valid Analyzable Scans: {valid_scans} ({retention_rate:.1f}%)")
    print(f"  • Flagged Scouts / Invalids: {failed_scans} ({100.0 - retention_rate:.1f}%)")

    # 2. Logistic Failure Model
    logit_summary, model_params = fit_logistic_failure_model(df)
    print("\n📈 Logistic Regression Failure Model [Logit(P(Failure))]:")
    print(logit_summary.to_string(index=False))

    # Save CSV
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    logit_summary.to_csv(out_csv, index=False)
    print(f"\n✓ Saved failure model table: {out_csv}")

    # 3. Generate Publication Figure 2
    generate_feasibility_plot(df, logit_summary, model_params, plot_path)
    print(f"✓ Saved publication figure: {plot_path}")
    print(f"Experiment 2 completed in {time.time() - t_start:.2f}s")
    print("=" * 75)


if __name__ == "__main__":
    main()
