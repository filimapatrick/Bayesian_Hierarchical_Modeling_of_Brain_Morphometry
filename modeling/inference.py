#!/usr/bin/env python3
"""
PyMC NUTS Inference Engine for Bayesian Brain Morphometry

Features:
- Full Hamiltonian Monte Carlo (HMC) with the No-U-Turn Sampler (NUTS) via PyMC 5.x
- Non-centered hierarchical parameterization for robust geometry and divergence-free sampling
- Heteroskedastic exponential slice noise dispersion: σ_i = σ_0 * exp(λ * (h_i - 1.0))
- Contrast adjustment parameter (δ)
- Site-level random intercepts with sum-to-zero constraint (γ_s)
- Clinical diagnostic random effects with reference control (β_d, β_CONTROL = 0)
- Modern convergence diagnostics via ArviZ: rank-normalized split R-hat, bulk ESS, tail ESS, 95% HDI
- Variance partitioning: ICC_site (Site-level variance fraction) vs ICC_disorder (Disorder share)
- Exports structured CSV summary tables and compressed posterior traces
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path
from typing import Dict, List, Tuple

# Ensure scipy.signal.gaussian compatibility for ArviZ under SciPy 1.13+
import scipy.signal
import scipy.signal.windows
if not hasattr(scipy.signal, "gaussian"):
    scipy.signal.gaussian = scipy.signal.windows.gaussian

import arviz as az
import numpy as np
import pandas as pd
import pymc as pm

PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))


def canonicalize_site(name: str) -> str:
    """Canonicalizes raw hospital institution names into unified canonical site IDs."""
    if not isinstance(name, str):
        return "UNKNOWN"
    upper = name.strip().upper()
    if "RSUTH" in upper:
        return "RSUTH"
    elif "UPTH" in upper:
        return "UPTH"
    elif "AKTH" in upper or "NKDC" in upper or "AMINU" in upper:
        return "AKTH_NKDC"
    elif "BMH" in upper or "BRAITHWAITE" in upper:
        return "BMH"
    elif "LIFE" in upper or "BRIDGE" in upper:
        return "LifeBridge"
    elif "IDC" in upper or "INTERCONTINENTAL" in upper:
        return "IDC"
    return name.strip()


def build_pymc_model(
    df: pd.DataFrame,
    target_metric: str = "evans_index",
    priors_config: Dict[str, float] = None,
) -> Tuple[pm.Model, List[str], List[str]]:
    """
    Constructs the generative PyMC hierarchical model.
    """
    if priors_config is None:
        if target_metric.lower() == "evans_index":
            priors_config = {"alpha_mu": 0.60, "alpha_sd": 0.10, "sigma_0_sd": 0.05}
        elif target_metric.lower() in ["bpf", "pef"]:
            priors_config = {"alpha_mu": 0.75, "alpha_sd": 0.15, "sigma_0_sd": 0.05}
        else:  # VBR
            priors_config = {"alpha_mu": 0.01, "alpha_sd": 0.01, "sigma_0_sd": 0.01}

    # Extract clean vectors
    y = df[target_metric].values.astype(float)
    h = df["slice_thickness_mm"].values.astype(float)
    c = df["contrast_enhanced"].astype(float).values

    # Encode diagnoses with CONTROL as reference (beta_CONTROL = 0)
    all_diagnoses = sorted(df["diagnosis"].unique())
    ref_diag = "CONTROL"
    non_ref_diags = [d for d in all_diagnoses if d != ref_diag]

    X_diag = np.zeros((len(df), len(non_ref_diags)), dtype=float)
    for i, d in enumerate(df["diagnosis"]):
        if d in non_ref_diags:
            X_diag[i, non_ref_diags.index(d)] = 1.0

    # Encode sites
    sites = sorted(df["site"].unique())
    site_idx = np.array([sites.index(s) for s in df["site"]], dtype=int)
    num_sites = len(sites)

    with pm.Model() as model:
        # 1. Global baseline & contrast shift
        alpha = pm.Normal("alpha", mu=priors_config["alpha_mu"], sigma=priors_config["alpha_sd"])
        delta = pm.Normal("delta", mu=0.0, sigma=0.05)

        # 2. Heteroskedastic exponential noise dispersion
        sigma_0 = pm.HalfNormal("sigma_0", sigma=priors_config["sigma_0_sd"])
        lam = pm.Normal("lambda", mu=0.0, sigma=0.15)
        sigma_i = pm.Deterministic("sigma_i", sigma_0 * pm.math.exp(lam * (h - 1.0)))

        # 3. Clinical disorder random effects (non-centered partial pooling)
        tau_beta = pm.HalfNormal("tau_beta", sigma=0.05)
        beta_raw = pm.Normal("beta_raw", mu=0.0, sigma=1.0, shape=len(non_ref_diags))
        beta_non_ref = pm.Deterministic("beta_non_ref", beta_raw * tau_beta)

        # 4. Site random intercepts (non-centered with sum-to-zero constraint)
        tau_gamma = pm.HalfNormal("tau_gamma", sigma=0.05)
        gamma_raw = pm.Normal("gamma_raw", mu=0.0, sigma=1.0, shape=num_sites)
        gamma_zero_mean = gamma_raw - pm.math.mean(gamma_raw)
        gamma = pm.Deterministic("gamma", gamma_zero_mean * tau_gamma)

        # 5. Linear predictor
        mu_i = alpha + delta * c + gamma[site_idx] + pm.math.dot(X_diag, beta_non_ref)

        # 6. Variance Partitioning / Intraclass Correlation Coefficients
        # Thin-slice (1.0 mm) idealized reference variance fraction:
        var_total_1mm = tau_gamma**2 + tau_beta**2 + sigma_0**2
        icc_site_1mm = pm.Deterministic("ICC_site_1mm", tau_gamma**2 / var_total_1mm)
        icc_disorder_1mm = pm.Deterministic("ICC_disorder_1mm", tau_beta**2 / var_total_1mm)

        # Clinical cohort average heteroskedastic residual variance:
        # sigma_i^2 = sigma_0^2 * exp(2 * lambda * (h - 1.0))
        mean_sigma2_cohort = pm.math.mean(sigma_i**2)
        var_total_clinical = tau_gamma**2 + tau_beta**2 + mean_sigma2_cohort
        icc_site_clinical = pm.Deterministic("ICC_site_clinical", tau_gamma**2 / var_total_clinical)
        icc_disorder_clinical = pm.Deterministic("ICC_disorder_clinical", tau_beta**2 / var_total_clinical)

        # Legacy alias for backward compatibility
        icc_site = pm.Deterministic("ICC_site", icc_site_1mm)
        icc_disorder = pm.Deterministic("ICC_disorder", icc_disorder_1mm)

        # 7. Likelihood
        pm.Normal("obs", mu=mu_i, sigma=sigma_i, observed=y)

    return model, non_ref_diags, sites


def fit_bayesian_model(
    df: pd.DataFrame,
    target_metric: str = "evans_index",
    draws: int = 2000,
    tune: int = 1000,
    chains: int = 4,
    random_seed: int = 42,
) -> Tuple[az.InferenceData, pd.DataFrame]:
    """
    Fits the Bayesian model via PyMC NUTS and generates structured posterior summary.
    """
    df_clean = df.copy()
    if "site" not in df_clean.columns:
        df_clean["site"] = df_clean["institution_name"].apply(canonicalize_site)

    # Filter out invalid or missing measurements
    valid_mask = (df_clean["is_volumetric_valid"] == 1.0) & df_clean[target_metric].notna()
    df_fit = df_clean[valid_mask].copy()

    print(f"Loaded {len(df_fit)} valid subjects across {df_fit['diagnosis'].nunique()} clinical groups "
          f"and {df_fit['site'].nunique()} scanner sites.")
    print(f"Disorder Levels: {sorted(df_fit['diagnosis'].unique())}")
    print(f"Site Levels:     {sorted(df_fit['site'].unique())}")

    model, non_ref_diags, sites = build_pymc_model(df_fit, target_metric=target_metric)

    print(f"\nSampling with PyMC NUTS (Chains: {chains}, Draws: {draws}, Tune: {tune}, Target Accept: 0.95)...")
    t0 = time.time()
    with model:
        idata = pm.sample(
            draws=draws,
            tune=tune,
            chains=chains,
            target_accept=0.95,
            random_seed=random_seed,
            return_inferencedata=True,
            progressbar=True,
        )
    sampling_time = time.time() - t0
    print(f"✓ All {chains} chains completed in {sampling_time:.1f}s ({sampling_time/chains:.1f}s/chain)")

    # Extract ArviZ summary
    core_vars = [
        "alpha", "delta", "sigma_0", "lambda", "tau_beta", "tau_gamma",
        "ICC_site", "ICC_disorder", "ICC_site_clinical", "ICC_disorder_clinical"
    ]
    az_core = az.summary(idata, var_names=core_vars, hdi_prob=0.95)
    az_beta = az.summary(idata, var_names=["beta_non_ref"], hdi_prob=0.95)
    az_gamma = az.summary(idata, var_names=["gamma"], hdi_prob=0.95)

    # Format parameter names
    records = []
    
    def add_row(param_name, row):
        records.append({
            "Parameter": param_name,
            "Mean": float(row["mean"]),
            "SD": float(row["sd"]),
            "HDI_2.5%": float(row["hdi_2.5%"]),
            "HDI_97.5%": float(row["hdi_97.5%"]),
            "R_hat": float(row["r_hat"]),
            "ESS_bulk": float(row["ess_bulk"]),
            "ESS_tail": float(row["ess_tail"]),
        })

    # Core parameters
    add_row("alpha (Global Baseline)", az_core.loc["alpha"])
    add_row("delta (Contrast Effect)", az_core.loc["delta"])
    add_row("sigma_0 (Baseline Noise at 1mm)", az_core.loc["sigma_0"])
    add_row("lambda (Slice Noise Scale)", az_core.loc["lambda"])
    add_row("tau_beta (Disorder SD)", az_core.loc["tau_beta"])
    add_row("tau_gamma (Site SD)", az_core.loc["tau_gamma"])

    # Reference diagnosis
    records.append({
        "Parameter": "beta_CONTROL (Ref=0)",
        "Mean": 0.0,
        "SD": 0.0,
        "HDI_2.5%": 0.0,
        "HDI_97.5%": 0.0,
        "R_hat": 1.000,
        "ESS_bulk": float(draws * chains),
        "ESS_tail": float(draws * chains),
    })

    # Non-reference diagnoses
    for idx, d in enumerate(non_ref_diags):
        key = f"beta_non_ref[{idx}]"
        if key in az_beta.index:
            add_row(f"beta_{d}", az_beta.loc[key])

    # Sites
    for idx, s in enumerate(sites):
        key = f"gamma[{idx}]"
        if key in az_gamma.index:
            add_row(f"gamma_{s}", az_gamma.loc[key])

    # Variance decomposition (Dual reporting: 1mm idealized reference vs clinical cohort average)
    add_row("ICC_disorder (Disorder Share - 1mm Ref)", az_core.loc["ICC_disorder"])
    add_row("ICC_site (Site Variance Fraction - 1mm Ref)", az_core.loc["ICC_site"])
    add_row("ICC_disorder_clinical (Cohort Average Disorder Share)", az_core.loc["ICC_disorder_clinical"])
    add_row("ICC_site_clinical (Cohort Average Site Variance Fraction)", az_core.loc["ICC_site_clinical"])

    summary_df = pd.DataFrame(records)
    return idata, summary_df


def print_summary_table(summary_df: pd.DataFrame, target_metric: str):
    """Prints a publication-ready formatted ASCII table."""
    print("\n" + "=" * 95)
    print(f"POSTERIOR ESTIMATES & CONVERGENCE DIAGNOSTICS ({target_metric.upper()}) [PyMC NUTS]")
    print("=" * 95)
    print(f"{'Parameter':<38} {'Mean':>8} {'SD':>8} {'95% HDI':>22} {'R-hat':>8} {'ESS Bulk':>10}")
    print("-" * 95)
    for _, row in summary_df.iterrows():
        p_name = row["Parameter"]
        mean_val = row["Mean"]
        sd_val = row["SD"]
        lo = row["HDI_2.5%"]
        hi = row["HDI_97.5%"]
        rhat = row["R_hat"]
        ess = row["ESS_bulk"]
        hdi_str = f"[{lo:8.4f}, {hi:8.4f}]"
        print(f"{p_name:<38} {mean_val:8.4f} {sd_val:8.4f} {hdi_str:>22} {rhat:8.3f} {ess:10.0f}")
    print("=" * 95)


def main():
    parser = argparse.ArgumentParser(description="Run PyMC NUTS Bayesian Inference for Brain Morphometry.")
    parser.add_argument("--features_csv", type=str, default="results/tables/macro_features.csv")
    parser.add_argument("--target_metric", type=str, default="evans_index", choices=["bpf", "pef", "vbr", "evans_index"])
    parser.add_argument("--draws", type=int, default=2000)
    parser.add_argument("--tune", type=int, default=1000)
    parser.add_argument("--chains", type=int, default=4)
    parser.add_argument("--out_summary", type=str, default=None)
    parser.add_argument("--out_traces", type=str, default=None)
    args = parser.parse_args()

    feat_path = Path(args.features_csv)
    if not feat_path.exists():
        sys.exit(f"Error: {feat_path} does not exist. Run feature extraction first.")

    raw_metric = args.target_metric.lower()
    df = pd.read_csv(feat_path)

    # Resolve metric column: map pef <-> bpf
    if raw_metric in ["pef", "bpf"]:
        metric_col = "pef" if "pef" in df.columns else "bpf"
        canonical_name = "pef"
    else:
        metric_col = raw_metric
        canonical_name = raw_metric

    out_summary = Path(args.out_summary) if args.out_summary else Path(f"results/tables/posterior_summary_{canonical_name}.csv")
    out_traces = Path(args.out_traces) if args.out_traces else Path(f"results/traces/mcmc_traces_{canonical_name}.npz")

    out_summary.parent.mkdir(parents=True, exist_ok=True)
    out_traces.parent.mkdir(parents=True, exist_ok=True)

    print("=" * 70)
    print(f"Hierarchical Bayesian NUTS Sampling: Target Metric = {canonical_name.upper()} (column: {metric_col})")
    print(f"Engine: PyMC v{pm.__version__} | ArviZ v{az.__version__}")
    print(f"Chains: {args.chains} | Draws: {args.draws} | Tuning: {args.tune}")
    print("=" * 70)

    idata, summary_df = fit_bayesian_model(
        df=df,
        target_metric=metric_col,
        draws=args.draws,
        tune=args.tune,
        chains=args.chains,
    )

    # Save canonical summary table
    summary_df.to_csv(out_summary, index=False)
    print(f"\n✓ Saved posterior summary table: {out_summary}")

    # If pef, also maintain posterior_summary_bpf.csv for backward compatibility
    if canonical_name == "pef":
        bpf_summary = Path("results/tables/posterior_summary_bpf.csv")
        summary_df.to_csv(bpf_summary, index=False)
        print(f"✓ Mirrored to legacy path: {bpf_summary}")

    # Save posterior traces
    posterior_dict = {}
    for var_name in idata.posterior.data_vars:
        posterior_dict[var_name] = idata.posterior[var_name].values
    np.savez_compressed(out_traces, **posterior_dict)
    print(f"✓ Saved compressed posterior traces: {out_traces}")

    if canonical_name == "pef":
        bpf_traces = Path("results/traces/mcmc_traces_bpf.npz")
        np.savez_compressed(bpf_traces, **posterior_dict)
        print(f"✓ Mirrored to legacy path: {bpf_traces}")

    print_summary_table(summary_df, canonical_name)


if __name__ == "__main__":
    main()
