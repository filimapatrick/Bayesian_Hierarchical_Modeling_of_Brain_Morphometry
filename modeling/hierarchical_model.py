import pymc as pm
import numpy as np

def build_model(y, site_idx, disorder_idx, X):
    """
    Hierarchical Bayesian model for brain morphometry
    """

    n_sites = len(np.unique(site_idx))
    n_disorders = len(np.unique(disorder_idx))

    with pm.Model() as model:

        # Hyperpriors
        mu_alpha = pm.Normal("mu_alpha", mu=0, sigma=10)
        sigma_alpha = pm.HalfNormal("sigma_alpha", sigma=5)

        # Site effects
        site_effect = pm.Normal("site_effect", mu=0, sigma=sigma_alpha, shape=n_sites)

        # Disorder effects
        disorder_effect = pm.Normal("disorder_effect", mu=0, sigma=5, shape=n_disorders)

        # Covariates
        beta = pm.Normal("beta", mu=0, sigma=5, shape=X.shape[1])

        # Linear model
        mu = (
            mu_alpha
            + site_effect[site_idx]
            + disorder_effect[disorder_idx]
            + pm.math.dot(X, beta)
        )

        sigma = pm.HalfNormal("sigma", sigma=5)

        # Likelihood
        y_obs = pm.Normal("y_obs", mu=mu, sigma=sigma, observed=y)

    return model