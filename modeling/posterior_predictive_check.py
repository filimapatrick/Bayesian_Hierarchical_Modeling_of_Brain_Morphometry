#!/usr/bin/env python3
"""
Posterior Predictive Checks (PPC) for Hierarchical Bayesian Brain Morphometry

Evaluates generative model fidelity and empirical coverage:
1. Simulates replicated data: y_rep ~ p(y | theta)
2. Compares observed vs. replicated distributions (mean, variance, tails)
3. Evaluates boundary behavior for bounded metrics: P(y_rep < 0) and P(y_rep > 1)
4. Renders publication Figure 5: Posterior Predictive Distributions vs Observed Data
5. Exports structured PPC summary table: results/tables/posterior_predictive_summary.csv
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

PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from modeling.inference import canonicalize_site


def run_ppc_for_metric(
    df_valid: pd.DataFrame,
    target_metric: str,
    trace_path: Path,
    num_samples: int = 2000,
    random_seed: int = 42,
) -> Tuple[np.ndarray, np.ndarray, Dict[str, float]]:
    """
    Computes posterior predictive replications for a given morphometric metric.
    """
    if not trace_path.exists():
        raise FileNotFoundError(f"Trace file not found: {trace_path}")

    traces = np.load(trace_path)
    
    # Map diagnoses
    all_diagnoses = sorted(df_valid["diagnosis"].unique())
    ref_diag = "CONTROL"
    non_ref_diags = [d for d in all_diagnoses if d != ref_diag]
    
    # Map sites
    sites = sorted(df_valid["site"].unique())
    site_map = {s: i for i, s in enumerate(sites)}
    
    y_obs = df_valid[target_metric].values.astype(float)
    h = df_valid["slice_thickness_mm"].values.astype(float)
    c = df_valid["contrast_enhanced"].astype(float).values
    
    # Flatten traces across chains if needed
    alpha = traces["alpha"].flatten()
    delta = traces["delta"].flatten()
    lam = traces["lambda"].flatten()
    sigma_0 = traces["sigma_0"].flatten()
    
    beta_non_ref = traces["beta_non_ref"]
    if beta_non_ref.ndim == 3:
        n_chains, n_draws, n_params = beta_non_ref.shape
        beta_non_ref = beta_non_ref.reshape(-1, n_params)
    
    gamma = traces["gamma"]
    if gamma.ndim == 3:
        n_chains, n_draws, n_sites = gamma.shape
        gamma = gamma.reshape(-1, n_sites)

    total_draws = len(alpha)
    rng = np.random.default_rng(random_seed)
    sampled_indices = rng.choice(total_draws, size=min(num_samples, total_draws), replace=False)
    
    n_obs = len(y_obs)
    y_rep = np.zeros((len(sampled_indices), n_obs))
    
    for s_idx, draw_idx in enumerate(sampled_indices):
        a_s = alpha[draw_idx]
        d_s = delta[draw_idx]
        l_s = lam[draw_idx]
        s0_s = sigma_0[draw_idx]
        b_s = beta_non_ref[draw_idx]
        g_s = gamma[draw_idx]
        
        # Heteroskedastic noise per subject
        sigma_i_s = s0_s * np.exp(l_s * (h - 1.0))
        
        # Linear predictor
        mu_s = a_s + d_s * c
        for i, row in df_valid.reset_index(drop=True).iterrows():
            diag = row["diagnosis"]
            if diag in non_ref_diags:
                mu_s[i] += b_s[non_ref_diags.index(diag)]
            s_name = row["site"]
            if s_name in site_map:
                mu_s[i] += g_s[site_map[s_name]]
                
        y_rep[s_idx, :] = rng.normal(mu_s, sigma_i_s)
        
    y_rep_flat = y_rep.flatten()
    
    stats_dict = {
        "Metric": target_metric.upper(),
        "Observed_Mean": float(np.mean(y_obs)),
        "Observed_SD": float(np.std(y_obs)),
        "Observed_Min": float(np.min(y_obs)),
        "Observed_Max": float(np.max(y_obs)),
        "Replicated_Mean": float(np.mean(y_rep_flat)),
        "Replicated_SD": float(np.std(y_rep_flat)),
        "Replicated_Min": float(np.min(y_rep_flat)),
        "Replicated_Max": float(np.max(y_rep_flat)),
        "Pct_Replications_Below_Zero": float((y_rep_flat < 0).mean() * 100.0),
        "Pct_Replications_Above_One": float((y_rep_flat > 1.0).mean() * 100.0),
    }
    
    return y_obs, y_rep, stats_dict


def render_ppc_figure(
    results_dict: Dict[str, Tuple[np.ndarray, np.ndarray, Dict[str, float]]],
    output_path: Path,
):
    """
    Renders publication Figure 5:
    Posterior Predictive Density Overlays vs Observed Empirical Distributions.
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)
    sns.set_theme(style="whitegrid", font_scale=1.05)
    
    fig, axes = plt.subplots(1, 3, figsize=(18, 5.5))
    metrics_info = [
        ("evans_index", "A. Evans' Index (EI)", axes[0], "#1f77b4", (0.1, 0.9)),
        ("pef", "B. Parenchymal Envelope Fraction (PEF)", axes[1], "#2ca02c", (0.35, 1.15)),
        ("vbr", "C. Ventricle-to-Brain Ratio (VBR)", axes[2], "#d62728", (-0.03, 0.08)),
    ]
    
    for metric_key, title, ax, color, xlims in metrics_info:
        if metric_key not in results_dict:
            continue
        y_obs, y_rep, stats = results_dict[metric_key]
        
        # Plot 50 individual posterior predictive density realizations
        num_plot_draws = min(50, len(y_rep))
        for i in range(num_plot_draws):
            sns.kdeplot(y_rep[i, :], ax=ax, color=color, alpha=0.08, linewidth=1.0)
            
        # Plot average posterior predictive density
        sns.kdeplot(y_rep.flatten(), ax=ax, color="#333333", linestyle="--", linewidth=2.0, label="PPC Mean Replication")
        
        # Plot observed empirical distribution
        sns.kdeplot(y_obs, ax=ax, color=color, linewidth=2.8, label="Observed Clinical Data", fill=True, alpha=0.25)
        
        ax.set_title(title, fontweight="bold", pad=12)
        ax.set_xlabel(f"Measured {metric_key.upper()} Value", fontweight="bold")
        ax.set_ylabel("Probability Density", fontweight="bold")
        ax.set_xlim(xlims)
        
        # Annotate boundary percentages
        below_0 = stats["Pct_Replications_Below_Zero"]
        above_1 = stats["Pct_Replications_Above_One"]
        info_text = (
            f"Observed: $\\mu={stats['Observed_Mean']:.3f}, \\sigma={stats['Observed_SD']:.3f}$\n"
            f"Replicated: $\\mu={stats['Replicated_Mean']:.3f}, \\sigma={stats['Replicated_SD']:.3f}$\n"
            f"$P(y^{{rep}} < 0) = {below_0:.2f}\\%$\n"
            f"$P(y^{{rep}} > 1) = {above_1:.2f}\\%$"
        )
        ax.text(
            0.05, 0.92, info_text,
            transform=ax.transAxes,
            fontsize=8.5,
            verticalalignment="top",
            bbox=dict(boxstyle="round,pad=0.5", facecolor="white", edgecolor="#cccccc", alpha=0.9),
        )
        ax.legend(loc="upper right", fontsize=8.5, framealpha=0.9)
        
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()


def main():
    parser = argparse.ArgumentParser(description="Run Posterior Predictive Checks (PPC).")
    parser.add_argument("--features_csv", type=str, default="results/tables/macro_features.csv")
    parser.add_argument("--traces_dir", type=str, default="results/traces")
    parser.add_argument("--output_csv", type=str, default="results/tables/posterior_predictive_summary.csv")
    parser.add_argument("--output_plot", type=str, default="results/figures/figure5_posterior_predictive_checks.png")
    args = parser.parse_args()

    feat_path = Path(args.features_csv)
    traces_dir = Path(args.traces_dir)
    out_csv = Path(args.output_csv)
    out_plot = Path(args.output_plot)

    if not feat_path.exists():
        sys.exit(f"Error: {feat_path} does not exist.")

    df = pd.read_csv(feat_path)
    df_valid = df[df["is_volumetric_valid"] == 1.0].copy()
    df_valid["site"] = df_valid["institution_name"].apply(canonicalize_site)

    print("=" * 80)
    print("POSTERIOR PREDICTIVE CHECKS (PPC) & EMPIRICAL COVERAGE EVALUATION")
    print(f"Cohort: {len(df_valid)} valid scans from {feat_path}")
    print("=" * 80)

    results = {}
    summary_records = []

    metrics = [
        ("evans_index", "evans_index", traces_dir / "mcmc_traces_evans_index.npz"),
        ("pef", "bpf" if "pef" not in df_valid.columns else "pef", traces_dir / "mcmc_traces_pef.npz"),
        ("vbr", "vbr", traces_dir / "mcmc_traces_vbr.npz"),
    ]

    for label, col, trace_p in metrics:
        if not trace_p.exists() and label == "pef":
            trace_p = traces_dir / "mcmc_traces_bpf.npz"
            
        print(f"\nEvaluating {label.upper()} (column: {col}) from {trace_p.name}...")
        y_obs, y_rep, stats = run_ppc_for_metric(df_valid, col, trace_p)
        results[label] = (y_obs, y_rep, stats)
        summary_records.append(stats)
        
        print(f"  • Observed:   mean = {stats['Observed_Mean']:.4f}, sd = {stats['Observed_SD']:.4f}, range = [{stats['Observed_Min']:.4f}, {stats['Observed_Max']:.4f}]")
        print(f"  • Replicated: mean = {stats['Replicated_Mean']:.4f}, sd = {stats['Replicated_SD']:.4f}, range = [{stats['Replicated_Min']:.4f}, {stats['Replicated_Max']:.4f}]")
        print(f"  • Mass < 0:   {stats['Pct_Replications_Below_Zero']:.4f}%")
        print(f"  • Mass > 1:   {stats['Pct_Replications_Above_One']:.4f}%")

    summary_df = pd.DataFrame(summary_records)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    summary_df.to_csv(out_csv, index=False)
    print(f"\n✓ Saved PPC summary table: {out_csv}")

    render_ppc_figure(results, out_plot)
    print(f"✓ Saved publication Figure 5: {out_plot}")
    print("=" * 80)


if __name__ == "__main__":
    main()
