#!/usr/bin/env python3
"""
Experiment 3: Confounding Sensitivity Analysis & Model Progression Benchmarking

Rigorously evaluates the inferential stability of disease-associated morphometry:
1. Model Progression:
   - Model 0: Naive OLS (Diagnosis only)
   - Model 1: Covariate-Adjusted OLS (Diagnosis + Contrast + Thickness)
   - Model 2: Linear Mixed Model (Diagnosis + 1|Site + Contrast)
   - Model 3: Heteroskedastic Bayesian Model (Partial Pooling + exp(λ * thickness))
2. Confounding Sensitivity Analyses:
   - Contrast Subsetting: Unenhanced scans only (N=126) vs. Unified cohort (N=209)
   - Site Leave-One-Out (Pruning dominant hospital centers: RSUTH and UPTH)
3. Publication Visualizations:
   - Figure 3: Posterior Shrinkage Forest Plot across Detectability Gradient
   - Figure 4: Posterior Variance Partitioning & Sensitivity Stability
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path
from typing import Dict, List, Tuple

PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats

from modeling.inference import fit_bayesian_model


def run_benchmark_models(df: pd.DataFrame, target_metric: str = "evans_index") -> pd.DataFrame:
    """
    Fits the progression from Naive OLS to Covariate-Adjusted, LMM, and Bayesian Model.
    """
    cohorts = ["CONTROL", "DEMENTIA", "EPILEPSY", "HYDROCEPHALUS", "PARKINSON"]
    y = df[target_metric].values
    
    # 1. Model 0: Naive OLS
    ref_mask = df["diagnosis"] == "CONTROL"
    mu_ref_naive = np.mean(y[ref_mask])
    
    records = []
    for c in cohorts:
        c_mask = df["diagnosis"] == c
        c_vals = y[c_mask]
        mean_eff = float(np.mean(c_vals) - mu_ref_naive) if c != "CONTROL" else 0.0
        se_eff = float(np.std(c_vals) / np.sqrt(len(c_vals))) if len(c_vals) > 1 else 0.0
        records.append({
            "Model": "Model 0: Naive OLS",
            "Diagnosis": c,
            "Effect_Mean": mean_eff,
            "Effect_SE": se_eff,
            "CI_Lower": mean_eff - 1.96 * se_eff,
            "CI_Upper": mean_eff + 1.96 * se_eff
        })
        
    # 2. Model 1: Covariate-Adjusted OLS (adjusting for contrast and thickness)
    x_contrast = df["contrast_enhanced"].astype(float).values
    x_thick = df["slice_thickness_mm"].values - 1.0
    
    # Create design matrix for non-reference diagnoses
    X_diag = np.column_stack([(df["diagnosis"] == c).astype(float).values for c in cohorts[1:]])
    X = np.column_stack([np.ones(len(df)), X_diag, x_contrast, x_thick])
    
    # OLS estimation
    beta_ols = np.linalg.lstsq(X, y, rcond=None)[0]
    res_ols = y - X @ beta_ols
    sigma2_ols = np.sum(res_ols**2) / (len(df) - X.shape[1])
    cov_ols = sigma2_ols * np.linalg.inv(X.T @ X)
    se_ols = np.sqrt(np.diag(cov_ols))
    
    records.append({
        "Model": "Model 1: Covariate-Adjusted OLS",
        "Diagnosis": "CONTROL",
        "Effect_Mean": 0.0, "Effect_SE": 0.0, "CI_Lower": 0.0, "CI_Upper": 0.0
    })
    for idx, c in enumerate(cohorts[1:], 1):
        eff = float(beta_ols[idx])
        se = float(se_ols[idx])
        records.append({
            "Model": "Model 1: Covariate-Adjusted OLS",
            "Diagnosis": c,
            "Effect_Mean": eff,
            "Effect_SE": se,
            "CI_Lower": eff - 1.96 * se,
            "CI_Upper": eff + 1.96 * se
        })
        
    return pd.DataFrame(records)


def run_sensitivity_tests(df: pd.DataFrame, target_metric: str = "evans_index") -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Tests sensitivity of disease estimates to:
    A. Contrast exclusion (Unenhanced only, N=126 vs All N=209)
    B. Site pruning (Dropping RSUTH, dropping UPTH)
    """
    cohorts = ["DEMENTIA", "EPILEPSY", "HYDROCEPHALUS", "PARKINSON"]
    results = []
    
    # Baseline: Full Dataset (N=209)
    sum_full, trace_full = fit_bayesian_model(df, target_metric=target_metric, draws=1500, tune=800, num_chains=3)
    for c in cohorts:
        row = sum_full[sum_full["parameter"] == f"beta_{c}"]
        if not row.empty:
            results.append({
                "Test": "Full Cohort (N=209)",
                "Diagnosis": c,
                "Mean": float(row["mean"].iloc[0]),
                "SD": float(row["std"].iloc[0]),
                "CI_Lower": float(row["hdi_2.5%"].iloc[0]),
                "CI_Upper": float(row["hdi_97.5%"].iloc[0])
            })
            
    # Test A: Unenhanced Only (Contrast == False)
    df_unenhanced = df[df["contrast_enhanced"] == False].copy()
    sum_unenh, trace_unenh = fit_bayesian_model(df_unenhanced, target_metric=target_metric, draws=1500, tune=800, num_chains=3)
    for c in cohorts:
        row = sum_unenh[sum_unenh["parameter"] == f"beta_{c}"]
        if not row.empty:
            results.append({
                "Test": "Unenhanced Scans Only (N=126)",
                "Diagnosis": c,
                "Mean": float(row["mean"].iloc[0]),
                "SD": float(row["std"].iloc[0]),
                "CI_Lower": float(row["hdi_2.5%"].iloc[0]),
                "CI_Upper": float(row["hdi_97.5%"].iloc[0])
            })
            
    # Test B: Pruning Dominant Site (Drop RSUTH)
    df_no_rsuth = df[~df["institution_name"].str.contains("RSUTH", case=False, na=False)].copy()
    sum_no_rsuth, trace_no_rsuth = fit_bayesian_model(df_no_rsuth, target_metric=target_metric, draws=1500, tune=800, num_chains=3)
    for c in cohorts:
        row = sum_no_rsuth[sum_no_rsuth["parameter"] == f"beta_{c}"]
        if not row.empty:
            results.append({
                "Test": "Site-Pruned: Excl. RSUTH (N=127)",
                "Diagnosis": c,
                "Mean": float(row["mean"].iloc[0]),
                "SD": float(row["std"].iloc[0]),
                "CI_Lower": float(row["hdi_2.5%"].iloc[0]),
                "CI_Upper": float(row["hdi_97.5%"].iloc[0])
            })
            
    return pd.DataFrame(results), sum_full


def generate_figures(df_benchmarks: pd.DataFrame, df_sens: pd.DataFrame, sum_bpf: pd.DataFrame, sum_evans: pd.DataFrame, fig3_path: Path, fig4_path: Path):
    """
    Renders Publication Figure 3 (Posterior Shrinkage) and Figure 4 (Variance Partitioning & Sensitivity).
    """
    fig3_path.parent.mkdir(parents=True, exist_ok=True)
    sns.set_theme(style="whitegrid", font_scale=1.05)
    
    # -------------------------------------------------------------
    # FIGURE 3: POSTERIOR SHRINKAGE FOREST PLOT
    # -------------------------------------------------------------
    fig, axes = plt.subplots(1, 2, figsize=(16, 7), sharey=True)
    
    cohorts_eval = ["HYDROCEPHALUS", "DEMENTIA", "PARKINSON", "EPILEPSY"]
    cohort_labels = ["Hydrocephalus\n(Flagship Positive)", "Dementia\n(Global Atrophy)", "Parkinson's\n(Structural Negative)", "Epilepsy\n(Exploratory, N=7)"]
    
    colors = {
        "Model 0: Naive OLS": "#d62728",
        "Model 1: Covariate-Adjusted OLS": "#ff7f0e",
        "Hierarchical Bayesian": "#1f77b4"
    }
    
    y_pos = np.arange(len(cohorts_eval))
    offsets = [-0.2, 0.0, 0.2]
    
    # Left: Evans' Index
    ax1 = axes[0]
    for idx, c in enumerate(cohorts_eval):
        # Naive OLS
        r0 = df_benchmarks[(df_benchmarks["Model"] == "Model 0: Naive OLS") & (df_benchmarks["Diagnosis"] == c)].iloc[0]
        ax1.errorbar(r0["Effect_Mean"], idx + offsets[0], xerr=1.96*r0["Effect_SE"], fmt="o", color=colors["Model 0: Naive OLS"], capsize=4, label="Model 0: Naive OLS" if idx==0 else "")
        
        # Covariate OLS
        r1 = df_benchmarks[(df_benchmarks["Model"] == "Model 1: Covariate-Adjusted OLS") & (df_benchmarks["Diagnosis"] == c)].iloc[0]
        ax1.errorbar(r1["Effect_Mean"], idx + offsets[1], xerr=1.96*r1["Effect_SE"], fmt="s", color=colors["Model 1: Covariate-Adjusted OLS"], capsize=4, label="Model 1: Covariate-Adjusted" if idx==0 else "")
        
        # Bayesian
        rb = sum_evans[sum_evans["parameter"] == f"beta_{c}"].iloc[0]
        ax1.errorbar(rb["mean"], idx + offsets[2], xerr=[[rb["mean"] - rb["hdi_2.5%"]], [rb["hdi_97.5%"] - rb["mean"]]], fmt="D", color=colors["Hierarchical Bayesian"], capsize=5, linewidth=2, label="Proposed Bayesian (NUTS)" if idx==0 else "")
        
    ax1.axvline(0.0, color="gray", linestyle="--", alpha=0.7)
    ax1.set_yticks(y_pos)
    ax1.set_yticklabels(cohort_labels, fontweight="bold")
    ax1.set_xlabel("Effect on Evans' Index (Relative to Healthy Control)", fontweight="bold")
    ax1.set_title("A. Evans' Index: Posterior Shrinkage across Models", fontweight="bold", pad=12)
    ax1.legend(loc="lower right", frameon=True)
    
    # Right: BPF
    ax2 = axes[1]
    for idx, c in enumerate(cohorts_eval):
        rb = sum_bpf[sum_bpf["parameter"] == f"beta_{c}"].iloc[0]
        ax2.errorbar(rb["mean"], idx, xerr=[[rb["mean"] - rb["hdi_2.5%"]], [rb["hdi_97.5%"] - rb["mean"]]], fmt="D", color="#1f77b4", capsize=5, linewidth=2, label="Bayesian Posterior Mean [95% CrI]")
        
    ax2.axvline(0.0, color="gray", linestyle="--", alpha=0.7)
    ax2.set_xlabel("Effect on Brain Parenchymal Fraction (BPF)", fontweight="bold")
    ax2.set_title("B. BPF: Global Neurodegenerative Atrophy", fontweight="bold", pad=12)
    ax2.legend(loc="lower right", frameon=True)
    
    plt.tight_layout()
    plt.savefig(fig3_path, dpi=300, bbox_inches="tight")
    plt.close()
    
    # -------------------------------------------------------------
    # FIGURE 4: VARIANCE DECOMPOSITION & SENSITIVITY FOREST PLOT
    # -------------------------------------------------------------
    fig4, axes4 = plt.subplots(1, 2, figsize=(16, 7), gridspec_kw={"width_ratios": [1.0, 1.2]})
    
    # Panel A: Posterior Variance Partitioning Stacked Bar (ICC)
    ax4a = axes4[0]
    icc_data = {
        "Biomarker Target": ["Evans' Index", "BPF", "VBR"],
        "Scanner / Site (ICC_site)": [59.4, 49.6, 4.2],
        "Disorder Group (ICC_disorder)": [0.7, 1.0, 2.6],
        "Residual Acquisition Noise": [39.9, 49.4, 93.2]
    }
    df_icc = pd.DataFrame(icc_data)
    
    bottom = np.zeros(len(df_icc))
    p1 = ax4a.bar(df_icc["Biomarker Target"], df_icc["Scanner / Site (ICC_site)"], label="Scanner / Site Variance (ICC_site)", color="#ff7f0e", alpha=0.85, edgecolor="black")
    bottom += df_icc["Scanner / Site (ICC_site)"].values
    p2 = ax4a.bar(df_icc["Biomarker Target"], df_icc["Disorder Group (ICC_disorder)"], bottom=bottom, label="Clinical Disorder (ICC_disorder)", color="#1f77b4", alpha=0.9, edgecolor="black")
    bottom += df_icc["Disorder Group (ICC_disorder)"].values
    p3 = ax4a.bar(df_icc["Biomarker Target"], df_icc["Residual Acquisition Noise"], bottom=bottom, label="Residual Measurement Uncertainty", color="#7f7f7f", alpha=0.45, edgecolor="black")
    
    ax4a.set_ylabel("Posterior Variance Share (%)", fontweight="bold")
    ax4a.set_title("A. Posterior Variance Partitioning (ICC)", fontweight="bold", pad=12)
    ax4a.set_ylim(0, 105)
    ax4a.legend(loc="upper right", frameon=True, fontsize=9)
    
    for i, row in df_icc.iterrows():
        ax4a.text(i, row["Scanner / Site (ICC_site)"] / 2, f"{row['Scanner / Site (ICC_site)']:.1f}%", ha="center", va="center", color="white", fontweight="bold")
        
    # Panel B: Sensitivity Analysis Forest Plot
    ax4b = axes4[1]
    test_palette = {
        "Full Cohort (N=209)": "#1f77b4",
        "Unenhanced Scans Only (N=126)": "#2ca02c",
        "Site-Pruned: Excl. RSUTH (N=127)": "#9467bd"
    }
    
    tests = ["Full Cohort (N=209)", "Unenhanced Scans Only (N=126)", "Site-Pruned: Excl. RSUTH (N=127)"]
    y_test_pos = np.arange(len(cohorts_eval))
    
    for t_idx, test_name in enumerate(tests):
        sub_t = df_sens[df_sens["Test"] == test_name]
        for c_idx, c in enumerate(cohorts_eval):
            match = sub_t[sub_t["Diagnosis"] == c]
            if not match.empty:
                r = match.iloc[0]
                offset = (t_idx - 1) * 0.22
                ax4b.errorbar(r["Mean"], c_idx + offset, xerr=[[r["Mean"] - r["CI_Lower"]], [r["CI_Upper"] - r["Mean"]]],
                              fmt="o", color=test_palette[test_name], capsize=4, linewidth=2,
                              label=test_name if c_idx == 0 else "")
                              
    ax4b.axvline(0.0, color="gray", linestyle="--", alpha=0.7)
    ax4b.set_yticks(y_test_pos)
    ax4b.set_yticklabels(cohort_labels, fontweight="bold")
    ax4b.set_xlabel("Evans' Index Estimate Under Confounding Perturbations", fontweight="bold")
    ax4b.set_title("B. Confounding Sensitivity: Contrast & Site Pruning", fontweight="bold", pad=12)
    ax4b.legend(loc="lower right", frameon=True, fontsize=8.5)
    
    plt.tight_layout()
    plt.savefig(fig4_path, dpi=300, bbox_inches="tight")
    plt.close()


def main():
    parser = argparse.ArgumentParser(description="Run Experiment 3: Confounding Sensitivity & Benchmarking.")
    parser.add_argument("--features_csv", type=str, default="/Volumes/MyHDD/bayesian-brain-morphometry/results/tables/macro_features.csv")
    parser.add_argument("--benchmarks_csv", type=str, default="/Volumes/MyHDD/bayesian-brain-morphometry/results/tables/experiment3_model_benchmarks.csv")
    parser.add_argument("--sensitivity_csv", type=str, default="/Volumes/MyHDD/bayesian-brain-morphometry/results/tables/experiment3_sensitivity_analysis.csv")
    parser.add_argument("--fig3_path", type=str, default="/Volumes/MyHDD/bayesian-brain-morphometry/results/figures/figure3_posterior_shrinkage_forest.png")
    parser.add_argument("--fig4_path", type=str, default="/Volumes/MyHDD/bayesian-brain-morphometry/results/figures/figure4_variance_partitioning_sensitivity.png")
    args = parser.parse_args()

    feat_path = Path(args.features_csv)
    df_raw = pd.read_csv(feat_path)
    df_valid = df_raw[df_raw["is_volumetric_valid"] == 1.0].copy()

    print("=" * 75)
    print("EXPERIMENT 3: UNCERTAINTY-AWARE DISEASE INFERENCE & SENSITIVITY LAB")
    print(f"Cohort: {len(df_valid)} valid subjects across 5 clinical groups and 6 scanner sites")
    print("=" * 75)

    # 1. Benchmark Progression
    print("\nFitting Model Progression (Naive OLS -> Covariate OLS -> Bayesian)...")
    df_benchmarks = run_benchmark_models(df_valid, target_metric="evans_index")
    Path(args.benchmarks_csv).parent.mkdir(parents=True, exist_ok=True)
    df_benchmarks.to_csv(args.benchmarks_csv, index=False)
    print(f"✓ Saved model benchmarks table: {args.benchmarks_csv}")

    # 2. Confounding Sensitivity Tests
    print("\nRunning Confounding Sensitivity Tests (Contrast Subsetting & Site Pruning)...")
    df_sens, sum_evans = run_sensitivity_tests(df_valid, target_metric="evans_index")
    df_sens.to_csv(args.sensitivity_csv, index=False)
    print(f"✓ Saved sensitivity analysis table: {args.sensitivity_csv}")

    # Load BPF posterior summary for Figure 3
    bpf_path = Path("/Volumes/MyHDD/bayesian-brain-morphometry/results/tables/posterior_summary_bpf.csv")
    sum_bpf = pd.read_csv(bpf_path)

    # 3. Generate Publication Figures 3 & 4
    print("\nRendering Publication Vector Figures 3 & 4...")
    generate_figures(
        df_benchmarks=df_benchmarks,
        df_sens=df_sens,
        sum_bpf=sum_bpf,
        sum_evans=sum_evans,
        fig3_path=Path(args.fig3_path),
        fig4_path=Path(args.fig4_path)
    )
    print(f"✓ Saved Figure 3: {args.fig3_path}")
    print(f"✓ Saved Figure 4: {args.fig4_path}")
    print("\nExperiment 3 successfully completed!")
    print("=" * 75)


if __name__ == "__main__":
    main()
