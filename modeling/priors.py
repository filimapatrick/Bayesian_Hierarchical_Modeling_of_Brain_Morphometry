#!/usr/bin/env python3
"""
Priors and Hyperprior Specifications for Bayesian Brain Morphometry

Defines weakly informative, regularizing prior distributions grounded in
biological plausibility and neuroanatomical scaling:
- Global intercept: alpha ~ Normal(mu_0, sigma_0^2)
- Disease group random effects: beta_k ~ Normal(0, tau_beta^2)
- Acquisition site random effects: gamma_j ~ Normal(0, tau_gamma^2)
- Hyperpriors on variance: tau_beta, tau_gamma ~ HalfNormal
- Contrast fixed effect: delta ~ Normal(0, 0.1)
- Heteroskedastic acquisition noise: sigma_i = sigma_0 * exp(lambda * thickness_i)
"""

from __future__ import annotations
from dataclasses import dataclass
import numpy as np


@dataclass(frozen=True)
class MorphometryPriors:
    # Baseline prior for target metric
    target_metric: str
    mu_alpha: float = 0.80      # Expected mean BPF for normal population
    sigma_alpha: float = 0.15   # Uncertainty on population mean
    
    # Hyperprior scale on disorder effects (between-group variability)
    tau_beta_scale: float = 0.20
    
    # Hyperprior scale on site effects (between-scanner variability)
    tau_gamma_scale: float = 0.10
    
    # Fixed effect prior scale for contrast enhancement
    delta_scale: float = 0.10
    
    # Baseline observation noise prior
    sigma_0_scale: float = 0.10
    
    # Slice thickness noise scaling coefficient
    lambda_mean: float = 0.0
    lambda_scale: float = 0.20

    @classmethod
    def get_priors_for_metric(cls, metric: str) -> MorphometryPriors:
        metric = metric.lower()
        if metric == "bpf":
            # Brain Parenchymal Fraction: range 0.70 - 0.95
            return cls(target_metric="bpf", mu_alpha=0.80, sigma_alpha=0.15, tau_beta_scale=0.15)
        elif metric == "vbr":
            # Ventricle-to-Brain Ratio: range 0.001 - 0.08
            return cls(target_metric="vbr", mu_alpha=0.010, sigma_alpha=0.02, tau_beta_scale=0.02)
        elif metric == "evans_index":
            # Evans' Index: range 0.40 - 0.75
            return cls(target_metric="evans_index", mu_alpha=0.60, sigma_alpha=0.15, tau_beta_scale=0.15)
        else:
            return cls(target_metric=metric, mu_alpha=0.0, sigma_alpha=1.0, tau_beta_scale=0.5)
