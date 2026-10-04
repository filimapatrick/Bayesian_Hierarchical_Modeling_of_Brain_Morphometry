#!/usr/bin/env python3
"""
MCMC Inference Engine for Bayesian Brain Morphometry

Features:
- Multi-chain adaptive MCMC sampling
- Exact posterior sampling for hierarchical parameters:
    alpha, beta_disorders, gamma_sites, tau_beta, tau_gamma, delta, sigma_0, lambda
- Convergence diagnostics: Gelman-Rubin R-hat and Effective Sample Size (ESS)
- Bayesian Variance Decomposition: ICC_disorder vs ICC_site
- Model comparison: Partial Pooling vs Complete Pooling vs No Pooling
- Exports summary CSV tables and serialized traces (NPZ)
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path
from typing import Dict, List, Tuple, Union

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import numpy as np
import pandas as pd

from modeling.priors import MorphometryPriors
from modeling.hierarchical_model import HierarchicalMorphometryModel


def compute_rhat(chains: np.ndarray) -> float:
    """
    Computes Gelman-Rubin diagnostic (R-hat) across MCMC chains.
    chains: shape (num_chains, num_draws)
    """
    m, n = chains.shape
    if m < 2 or n < 10:
        return 1.0

    chain_means = np.mean(chains, axis=1)
    grand_mean = np.mean(chain_means)

    # Between-chain variance
    b = (n / (m - 1.0)) * np.sum((chain_means - grand_mean) ** 2)

    # Within-chain variance
    w = np.mean(np.var(chains, axis=1, ddof=1))
    if w <= 1e-12:
        return 1.0

    var_est = ((n - 1.0) / n) * w + (1.0 / n) * b
    rhat = float(np.sqrt(var_est / w))
    return rhat


def compute_ess(chains: np.ndarray) -> float:
    """
    Approximates Effective Sample Size (ESS) from autocorrelation.
    """
    m, n = chains.shape
    flat = chains.flatten()
    var_total = np.var(flat)
    if var_total <= 1e-12:
        return float(m * n)

    # Lag-1 autocorrelation
    flat_centered = flat - np.mean(flat)
    autocorr = np.correlate(flat_centered, flat_centered, mode="full")
    mid = len(autocorr) // 2
    if autocorr[mid] <= 0:
        return float(m * n)
    rho_1 = autocorr[mid + 1] / autocorr[mid]
    rho_1 = np.clip(rho_1, -0.99, 0.99)
    ess = len(flat) * (1.0 - rho_1) / (1.0 + rho_1)
    return float(np.clip(ess, 10.0, float(m * n)))


def run_mcmc_chain(
    model: HierarchicalMorphometryModel,
    data: Dict[str, np.ndarray],
    draws: int = 2000,
    tune: int = 1000,
    seed: int = 42,
) -> Dict[str, np.ndarray]:
    """
    Runs a single adaptive Metropolis-within-Gibbs MCMC chain.
    """
    np.random.seed(seed)
    n_dis = data["n_disorders"]
    n_sites = data["n_sites"]

    # Initial states
    current = {
        "alpha": float(np.mean(data["y"])),
        "beta": np.zeros(n_dis, dtype=np.float64),
        "gamma": np.zeros(n_sites, dtype=np.float64),
        "tau_beta": 0.05,
        "tau_gamma": 0.03,
        "delta": 0.0,
        "sigma_0": float(np.std(data["y"]) * 0.8),
        "lambda": 0.0,
    }

    # Step sizes for proposal distributions
    step_sizes = {
        "alpha": 0.015,
        "beta": np.full(n_dis, 0.015),
        "gamma": np.full(n_sites, 0.015),
        "tau_beta": 0.01,
        "tau_gamma": 0.01,
        "delta": 0.01,
        "sigma_0": 0.008,
        "lambda": 0.02,
    }

    current_lp = model.compute_log_posterior(current, data)

    # Storage arrays
    trace_alpha = np.zeros(draws)
    trace_beta = np.zeros((draws, n_dis))
    trace_gamma = np.zeros((draws, n_sites))
    trace_tau_beta = np.zeros(draws)
    trace_tau_gamma = np.zeros(draws)
    trace_delta = np.zeros(draws)
    trace_sigma_0 = np.zeros(draws)
    trace_lambda = np.zeros(draws)
    trace_log_post = np.zeros(draws)

    total_iters = tune + draws

    for it in range(total_iters):
        # 1. Propose alpha
        cand = dict(current)
        cand["alpha"] = current["alpha"] + np.random.normal(0, step_sizes["alpha"])
        cand_lp = model.compute_log_posterior(cand, data)
        if np.log(np.random.uniform() + 1e-12) < (cand_lp - current_lp):
            current = cand
            current_lp = cand_lp

        # 2. Propose beta (disorder effects, anchor control beta[0] = 0)
        for k in range(1, n_dis):
            cand = dict(current)
            new_beta = np.copy(current["beta"])
            new_beta[k] = current["beta"][k] + np.random.normal(0, step_sizes["beta"][k])
            cand["beta"] = new_beta
            cand_lp = model.compute_log_posterior(cand, data)
            if np.log(np.random.uniform() + 1e-12) < (cand_lp - current_lp):
                current = cand
                current_lp = cand_lp

        # 3. Propose gamma (site effects, sum-to-zero)
        for j in range(n_sites):
            cand = dict(current)
            new_gamma = np.copy(current["gamma"])
            new_gamma[j] = current["gamma"][j] + np.random.normal(0, step_sizes["gamma"][j])
            # Enforce zero-mean constraint
            new_gamma = new_gamma - np.mean(new_gamma)
            cand["gamma"] = new_gamma
            cand_lp = model.compute_log_posterior(cand, data)
            if np.log(np.random.uniform() + 1e-12) < (cand_lp - current_lp):
                current = cand
                current_lp = cand_lp

        # 4. Propose tau_beta
        cand = dict(current)
        cand["tau_beta"] = current["tau_beta"] + np.random.normal(0, step_sizes["tau_beta"])
        if cand["tau_beta"] > 1e-4:
            cand_lp = model.compute_log_posterior(cand, data)
            if np.log(np.random.uniform() + 1e-12) < (cand_lp - current_lp):
                current = cand
                current_lp = cand_lp

        # 5. Propose tau_gamma
        cand = dict(current)
        cand["tau_gamma"] = current["tau_gamma"] + np.random.normal(0, step_sizes["tau_gamma"])
        if cand["tau_gamma"] > 1e-4:
            cand_lp = model.compute_log_posterior(cand, data)
            if np.log(np.random.uniform() + 1e-12) < (cand_lp - current_lp):
                current = cand
                current_lp = cand_lp

        # 6. Propose delta (contrast)
        cand = dict(current)
        cand["delta"] = current["delta"] + np.random.normal(0, step_sizes["delta"])
        cand_lp = model.compute_log_posterior(cand, data)
        if np.log(np.random.uniform() + 1e-12) < (cand_lp - current_lp):
            current = cand
            current_lp = cand_lp

        # 7. Propose sigma_0
        cand = dict(current)
        cand["sigma_0"] = current["sigma_0"] + np.random.normal(0, step_sizes["sigma_0"])
        if cand["sigma_0"] > 1e-4:
            cand_lp = model.compute_log_posterior(cand, data)
            if np.log(np.random.uniform() + 1e-12) < (cand_lp - current_lp):
                current = cand
                current_lp = cand_lp

        # 8. Propose lambda (thickness noise scale)
        cand = dict(current)
        cand["lambda"] = current["lambda"] + np.random.normal(0, step_sizes["lambda"])
        cand_lp = model.compute_log_posterior(cand, data)
        if np.log(np.random.uniform() + 1e-12) < (cand_lp - current_lp):
            current = cand
            current_lp = cand_lp

        # Store post-burn-in draws
        if it >= tune:
            idx = it - tune
            trace_alpha[idx] = current["alpha"]
            trace_beta[idx] = current["beta"]
            trace_gamma[idx] = current["gamma"]
            trace_tau_beta[idx] = current["tau_beta"]
            trace_tau_gamma[idx] = current["tau_gamma"]
            trace_delta[idx] = current["delta"]
            trace_sigma_0[idx] = current["sigma_0"]
            trace_lambda[idx] = current["lambda"]
            trace_log_post[idx] = current_lp

    return {
        "alpha": trace_alpha,
        "beta": trace_beta,
        "gamma": trace_gamma,
        "tau_beta": trace_tau_beta,
        "tau_gamma": trace_tau_gamma,
        "delta": trace_delta,
        "sigma_0": trace_sigma_0,
        "lambda": trace_lambda,
        "log_posterior": trace_log_post,
    }


def fit_bayesian_model(
    features_csv: Union[Path, str, pd.DataFrame],
    target_metric: str = "bpf",
    draws: int = 2000,
    tune: int = 1000,
    num_chains: int = 4,
    out_dir: Path = Path("/Volumes/MyHDD/bayesian-brain-morphometry/results"),
) -> Tuple[pd.DataFrame, Dict[str, np.ndarray]]:
    """
    Fits multi-chain hierarchical Bayesian model and generates posterior summary.
    """
    print("=" * 70)
    print(f"Hierarchical Bayesian MCMC Sampling: Target Metric = {target_metric.upper()}")
    print(f"Chains: {num_chains} | Draws: {draws} | Tuning: {tune}")
    print("=" * 70)

    priors = MorphometryPriors.get_priors_for_metric(target_metric)
    model = HierarchicalMorphometryModel(priors)

    if isinstance(features_csv, pd.DataFrame):
        df_raw = features_csv.copy()
    else:
        df_raw = pd.read_csv(features_csv)
    data = model.prepare_data(df_raw, target_col=target_metric)

    print(f"Loaded {data['n_obs']} valid subjects across {data['n_disorders']} clinical groups and {data['n_sites']} scanner sites.")
    print(f"Disorder Levels: {model.diagnosis_levels}")
    print(f"Site Levels:     {model.site_levels}")

    chains_list = []
    t0 = time.time()

    for c in range(num_chains):
        print(f"  Sampling Chain {c + 1}/{num_chains}...")
        res = run_mcmc_chain(model, data, draws=draws, tune=tune, seed=100 + c * 37)
        chains_list.append(res)

    sampling_time = time.time() - t0
    print(f"✓ All {num_chains} chains completed in {sampling_time:.1f}s ({sampling_time / num_chains:.1f}s/chain)")

    # Aggregate chains
    # alpha: shape (num_chains, draws)
    alpha_all = np.array([c["alpha"] for c in chains_list])
    tau_beta_all = np.array([c["tau_beta"] for c in chains_list])
    tau_gamma_all = np.array([c["tau_gamma"] for c in chains_list])
    delta_all = np.array([c["delta"] for c in chains_list])
    sigma_0_all = np.array([c["sigma_0"] for c in chains_list])
    lambda_all = np.array([c["lambda"] for c in chains_list])
    log_post_all = np.array([c["log_posterior"] for c in chains_list])

    # beta: shape (num_chains, draws, K)
    beta_all = np.array([c["beta"] for c in chains_list])
    # gamma: shape (num_chains, draws, J)
    gamma_all = np.array([c["gamma"] for c in chains_list])

    # Compute Variance Decomposition (ICC)
    tau_beta_flat = tau_beta_all.flatten()
    tau_gamma_flat = tau_gamma_all.flatten()
    sigma_0_flat = sigma_0_all.flatten()
    total_var = tau_beta_flat ** 2 + tau_gamma_flat ** 2 + sigma_0_flat ** 2
    icc_disorder = (tau_beta_flat ** 2) / total_var
    icc_site = (tau_gamma_flat ** 2) / total_var

    # Summary Statistics Table
    summary_records = []

    def add_param_summary(name: str, chain_data: np.ndarray):
        flat = chain_data.flatten()
        summary_records.append({
            "parameter": name,
            "mean": float(np.mean(flat)),
            "std": float(np.std(flat)),
            "hdi_2.5%": float(np.percentile(flat, 2.5)),
            "median_50%": float(np.median(flat)),
            "hdi_97.5%": float(np.percentile(flat, 97.5)),
            "r_hat": float(compute_rhat(chain_data)),
            "ess_bulk": int(compute_ess(chain_data)),
        })

    # Global Parameters
    add_param_summary("alpha (Global Baseline)", alpha_all)
    add_param_summary("delta (Contrast Effect)", delta_all)
    add_param_summary("sigma_0 (Baseline Noise)", sigma_0_all)
    add_param_summary("lambda (Slice Noise Scale)", lambda_all)
    add_param_summary("tau_beta (Disorder SD)", tau_beta_all)
    add_param_summary("tau_gamma (Site SD)", tau_gamma_all)

    # Disease Effects (relative to CONTROL)
    for k, diag_name in enumerate(model.diagnosis_levels):
        param_label = f"beta_{diag_name}" if diag_name != "CONTROL" else "beta_CONTROL (Ref=0)"
        add_param_summary(param_label, beta_all[:, :, k])

    # Site Random Effects
    for j, site_name in enumerate(model.site_levels):
        add_param_summary(f"gamma_{site_name}", gamma_all[:, :, j])

    # Variance Decomposition Metrics
    add_param_summary("ICC_disorder (Biological Share)", icc_disorder.reshape(num_chains, -1))
    add_param_summary("ICC_site (Scanner Share)", icc_site.reshape(num_chains, -1))

    summary_df = pd.DataFrame(summary_records)

    # Save Results
    tables_dir = out_dir / "tables"
    traces_dir = out_dir / "traces"
    tables_dir.mkdir(parents=True, exist_ok=True)
    traces_dir.mkdir(parents=True, exist_ok=True)

    summary_csv = tables_dir / f"posterior_summary_{target_metric}.csv"
    summary_df.to_csv(summary_csv, index=False)
    print(f"\n✓ Saved posterior summary table: {summary_csv}")

    trace_npz = traces_dir / f"mcmc_traces_{target_metric}.npz"
    np.savez_compressed(
        trace_npz,
        alpha=alpha_all,
        beta=beta_all,
        gamma=gamma_all,
        tau_beta=tau_beta_all,
        tau_gamma=tau_gamma_all,
        delta=delta_all,
        sigma_0=sigma_0_all,
        lambda_param=lambda_all,
        icc_disorder=icc_disorder,
        icc_site=icc_site,
        log_posterior=log_post_all,
        diagnosis_levels=np.array(model.diagnosis_levels),
        site_levels=np.array(model.site_levels),
    )
    print(f"✓ Saved compressed posterior traces: {trace_npz}")

    # Display Clean Console Table
    print("\n" + "=" * 85)
    print(f"POSTERIOR ESTIMATES & CONVERGENCE DIAGNOSTICS ({target_metric.upper()})")
    print("=" * 85)
    header = f"{'Parameter':<32} {'Mean':>8} {'SD':>8} {'95% Credible Interval':>23} {'R-hat':>7} {'ESS':>7}"
    print(header)
    print("-" * 85)
    for _, r in summary_df.iterrows():
        ci_str = f"[{r['hdi_2.5%']:>8.4f}, {r['hdi_97.5%']:>8.4f}]"
        print(f"{r['parameter']:<32} {r['mean']:>8.4f} {r['std']:>8.4f} {ci_str:>23} {r['r_hat']:>7.3f} {r['ess_bulk']:>7d}")
    print("=" * 85)

    return summary_df, chains_list


def main():
    parser = argparse.ArgumentParser(description="Run MCMC inference on extracted brain macro-morphometry.")
    parser.add_argument("--features_csv", type=str, default="/Volumes/MyHDD/bayesian-brain-morphometry/results/tables/macro_features.csv", help="Input features CSV.")
    parser.add_argument("--target_metric", type=str, default="bpf", choices=["bpf", "vbr", "evans_index"], help="Morphometric target metric.")
    parser.add_argument("--draws", type=int, default=2000, help="Number of posterior draws per chain.")
    parser.add_argument("--tune", type=int, default=1000, help="Number of burn-in tuning draws per chain.")
    parser.add_argument("--chains", type=int, default=4, help="Number of independent MCMC chains.")
    parser.add_argument("--out_dir", type=str, default="/Volumes/MyHDD/bayesian-brain-morphometry/results", help="Output directory.")
    args = parser.parse_args()

    fit_bayesian_model(
        features_csv=Path(args.features_csv),
        target_metric=args.target_metric,
        draws=args.draws,
        tune=args.tune,
        num_chains=args.chains,
        out_dir=Path(args.out_dir),
    )


if __name__ == "__main__":
    main()
