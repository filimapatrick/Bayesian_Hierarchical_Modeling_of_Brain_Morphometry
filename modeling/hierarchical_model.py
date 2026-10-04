#!/usr/bin/env python3
"""
Hierarchical Bayesian Model for Macro-Brain Morphometry

Mathematical Formulation:
y_i ~ Normal(mu_i, sigma_i^2)

mu_i = alpha + beta_{disorder[i]} + gamma_{site[i]} + delta * contrast_i
sigma_i = sigma_0 * exp(lambda * centered_slice_thickness_i)

Priors:
alpha ~ Normal(mu_alpha, sigma_alpha^2)
beta_k ~ Normal(0, tau_beta^2) with sum(beta) = 0 or relative to CONTROL
gamma_j ~ Normal(0, tau_gamma^2) with sum(gamma) = 0
tau_beta ~ HalfNormal(scale_beta)
tau_gamma ~ HalfNormal(scale_gamma)
delta ~ Normal(0, scale_delta^2)
sigma_0 ~ HalfNormal(scale_sigma_0)
lambda ~ Normal(0, scale_lambda^2)
"""

from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import numpy as np
import pandas as pd
from typing import Dict, Tuple, List, Optional
from modeling.priors import MorphometryPriors


class HierarchicalMorphometryModel:
    def __init__(self, priors: MorphometryPriors):
        self.priors = priors
        self.diagnosis_levels: List[str] = []
        self.site_levels: List[str] = []
        
    def prepare_data(self, df: pd.DataFrame, target_col: str = "bpf") -> Dict[str, np.ndarray]:
        # Filter valid records
        clean = df[df["is_volumetric_valid"] == 1.0].dropna(subset=[target_col]).copy()
        
        # Standardize site names (merge slight variations of RSUTH and Life Bridge)
        site_map = {
            "RSUTH Port Harcourt": "RSUTH",
            "RSUTH": "RSUTH",
            "UPTH": "UPTH",
            "INTERCONTINENTAL DIAG. CENTER": "IDC",
            "Aminu Kano Teaching Hospital / NKDC Kano": "AKTH_NKDC",
            "Braithwaite Memorial Hospital": "BMH",
            "Life Bridge": "LifeBridge",
            "LIFEBRIDGE MEDICAL DIAGNOSTICS LTD": "LifeBridge"
        }
        clean["clean_site"] = clean["institution_name"].map(lambda x: site_map.get(str(x).strip(), "Other"))
        
        # Categorical indexing
        # Make CONTROL the baseline index 0
        diag_list = sorted(clean["diagnosis"].unique())
        if "CONTROL" in diag_list:
            diag_list.remove("CONTROL")
            diag_list = ["CONTROL"] + diag_list
        self.diagnosis_levels = diag_list
        
        self.site_levels = sorted(clean["clean_site"].unique())
        
        clean["diag_idx"] = clean["diagnosis"].map(lambda x: self.diagnosis_levels.index(x))
        clean["site_idx"] = clean["clean_site"].map(lambda x: self.site_levels.index(x))
        clean["contrast_binary"] = clean["contrast_enhanced"].map(lambda x: 1.0 if str(x).lower() == "true" else 0.0)
        
        # Continuous covariates
        slice_thick = pd.to_numeric(clean["slice_thickness_mm"], errors="coerce").fillna(5.0).values
        mean_thick = np.mean(slice_thick)
        centered_thick = slice_thick - mean_thick
        
        y = clean[target_col].values.astype(np.float64)
        
        return {
            "y": y,
            "diag_idx": clean["diag_idx"].values.astype(np.int64),
            "site_idx": clean["site_idx"].values.astype(np.int64),
            "contrast": clean["contrast_binary"].values.astype(np.float64),
            "thickness": centered_thick.astype(np.float64),
            "raw_thickness": slice_thick,
            "n_obs": len(y),
            "n_disorders": len(self.diagnosis_levels),
            "n_sites": len(self.site_levels),
            "mean_thickness": mean_thick,
            "participant_ids": clean["participant_id"].values
        }

    def compute_log_posterior(self, params: Dict[str, np.ndarray], data: Dict[str, np.ndarray]) -> float:
        """
        Evaluates exact unnormalized log posterior density: log_prior + log_likelihood
        """
        alpha = params["alpha"]
        beta = params["beta"]            # shape (K,)
        gamma = params["gamma"]          # shape (J,)
        tau_beta = params["tau_beta"]
        tau_gamma = params["tau_gamma"]
        delta = params["delta"]
        sigma_0 = params["sigma_0"]
        lam = params["lambda"]
        
        # Positivity checks
        if tau_beta <= 1e-6 or tau_gamma <= 1e-6 or sigma_0 <= 1e-6:
            return -np.inf
            
        p = self.priors
        
        # 1. Log Priors
        lp = 0.0
        # alpha ~ Normal(mu_alpha, sigma_alpha^2)
        lp += -0.5 * ((alpha - p.mu_alpha) / p.sigma_alpha) ** 2
        
        # beta_k ~ Normal(0, tau_beta^2)
        lp += -0.5 * np.sum((beta / tau_beta) ** 2) - len(beta) * np.log(tau_beta)
        # tau_beta ~ HalfNormal(scale_beta)
        lp += -0.5 * (tau_beta / p.tau_beta_scale) ** 2
        
        # gamma_j ~ Normal(0, tau_gamma^2)
        lp += -0.5 * np.sum((gamma / tau_gamma) ** 2) - len(gamma) * np.log(tau_gamma)
        # tau_gamma ~ HalfNormal(scale_gamma)
        lp += -0.5 * (tau_gamma / p.tau_gamma_scale) ** 2
        
        # delta ~ Normal(0, delta_scale^2)
        lp += -0.5 * (delta / p.delta_scale) ** 2
        
        # sigma_0 ~ HalfNormal(sigma_0_scale)
        lp += -0.5 * (sigma_0 / p.sigma_0_scale) ** 2
        
        # lambda ~ Normal(lambda_mean, lambda_scale^2)
        lp += -0.5 * ((lam - p.lambda_mean) / p.lambda_scale) ** 2
        
        # 2. Log Likelihood
        diag_idx = data["diag_idx"]
        site_idx = data["site_idx"]
        contrast = data["contrast"]
        thickness = data["thickness"]
        y = data["y"]
        
        # Linear predictor mu_i
        # beta is anchored to CONTROL: beta[0] = 0.0 for identifiability
        mu = alpha + beta[diag_idx] + gamma[site_idx] + delta * contrast
        
        # Heteroskedastic noise model
        sigma_i = sigma_0 * np.exp(np.clip(lam * thickness, -3.0, 3.0))
        
        # Sum of Gaussian log likelihoods
        ll = -0.5 * np.sum(((y - mu) / sigma_i) ** 2) - np.sum(np.log(sigma_i))
        
        return float(lp + ll)