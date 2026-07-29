# bayesianOU

Bayesian Nonlinear Ornstein-Uhlenbeck Models with Stochastic Volatility

## Overview

`bayesianOU` fits a Bayesian nonlinear Ornstein-Uhlenbeck (OU) / nonlinear
error-correction model with cubic drift, stochastic volatility (SV), and
Student-t innovations. Sector-specific parameters use a **non-centered
hierarchical (partial-pooling)** parameterization, and sampling is done with
Stan (`reduce_sum` for within-chain parallelism). Model comparison uses
PSIS-LOO (Vehtari, Gelman, Gabry 2017).

## Installation

```r
# install.packages("remotes")
remotes::install_github("IsadoreNabi/bayesianOU")

# Stan backend (cmdstanr recommended):
install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
cmdstanr::install_cmdstan()
# Or rstan:
install.packages("rstan")
```

## Quick Start

```r
library(bayesianOU)

Y   <- as.matrix(your_prices_data)             # market prices by sector
X   <- as.matrix(your_production_prices_data)  # prices of production
TMG <- your_tmg_series                         # aggregate profit rate
COM <- as.matrix(your_com_data)                # organic composition of capital
K   <- as.matrix(your_capital_data)            # total capital by sector

results <- fit_ou_nonlinear_tmg(
  results_robust = list(),
  Y = Y, X = X, TMG = TMG, COM = COM, CAPITAL_TOTAL = K,
  fit_window = "train",   # honest out-of-sample evaluation
  chains = 4, iter = 4000, warmup = 2000,
  verbose = TRUE
)

validate_ou_fit(results)            # MCMC + LOO (incl. Pareto-k) + OOS
kappa_stability_evidence(results)   # dynamic mean-reversion evidence
plot_beta_tmg(results)
plot_drift_curves(results)
```

## Model Specification

For each sector `s`, on training-standardized series, the one-step (Euler,
dt = 1) increment is

```
dY[t,s] = kappa_s * (theta_s - Y[t-1,s] + a3_s * (Y[t-1,s] - theta_s)^3)   # cubic OU drift
        + (beta0_s + beta1 * zTMG[t]) * X[t-1,s]                           # TMG-modulated pass-through
        + gamma * COM_std[t-1,s]                                            # organic composition in mean
        + eps[t,s],     eps ~ Student-t(nu, 0, sigma[t,s])
log sigma^2[t,s] = h[t,s],   h ~ stationary AR(1)  (alpha_s, rho_s, sigma_eta_s)
```

- `kappa_s` mean-reversion speed, `theta_s` equilibrium level, `a3_s < 0` cubic
  term. The implied long-run equilibrium is `Y* = theta_s + (beta/kappa_s) X`,
  i.e. a nonlinear error-correction toward a linear function of the prices of
  production.
- `theta_s, kappa_s, a3_s, beta0_s` are built hierarchically as
  `intercept + slope * COM_s + sd * z`, `z ~ N(0,1)` (non-centered).
- `nu > 2` (finite variance). `sigma[t,s]` is the Student-t **scale**; the
  conditional standard deviation is `sigma * sqrt(nu/(nu-2))`.

## Methodology and validity (read before reporting results)

**Training/test and "out-of-sample".** A Bayesian fit on all data is perfectly
valid for *inference* and for *PSIS-LOO* (which approximates leave-one-out
predictive performance from the full-data posterior — this is exactly what
`rstanarm` does, and it is not leakage). Leakage only arises if one labels a
block "out-of-sample" while it still entered the likelihood. This package keeps
the two designs coherent via `fit_window`:

- `fit_window = "train"` (default): the likelihood and `log_lik` are summed
  **only** over the training window `2:T_train`; the test block is genuinely
  held out, so `evaluate_oos` is a real out-of-sample evaluation.
- `fit_window = "full"`: fit on all observations (full-information). PSIS-LOO is
  then valid over all points, but `evaluate_oos` becomes in-sample.

**Priors are configurable and neutral by default.** The key hypothesis prior is
neutral: `beta1 ~ Normal(0, 0.5)` (no sign baked in). The Student-t df uses a
weakly-informative `nu_tilde ~ Gamma(2, 0.1)` (prior mean `nu ~ 22`), and SV
persistence `rho_s ~ Normal(0.7, 0.2)`. Override any of these through the
`priors` argument and run a sensitivity analysis.

**Leakage in preprocessing is avoided by default.** `use_train_loadings = TRUE`
computes the common-factor loadings from the training window only and then
projects the full series, so the orthogonalized TMG regressor does not see the
future.

**Stability assumption.** `a3_s < 0` is enforced, so mean reversion strengthens
with the deviation. This is an *assumption* (global stabilization), not an
estimated result; it precludes detecting locally expansive / self-amplifying
regimes.

**Half-life.** With the Euler dt = 1 discretization, the discrete persistence of
the linear part is `1 - kappa`; interpret half-lives as
`log(0.5)/log(1 - kappa)`, not `log(2)/kappa` (these coincide only as
`kappa -> 0`).

**PSIS-LOO caveat.** The model has one latent SV state per observation, so plain
PSIS-LOO tends to be optimistic and to produce high Pareto-k. `validate_ou_fit`
surfaces a Pareto-k summary and warns; for genuine forecasting assessment prefer
leave-future-out (LFO-CV).

## Validation layers

- **Internal:** R-hat, ESS, divergences, Pareto-k (`validate_ou_fit`).
- **Synthetic recovery:** simulate from the generative process with known
  parameters, fit, and check coverage (see `tests/testthat/test-recovery.R`).
- **Organic / external:** held-out metrics (`fit_window = "train"`) and PSIS-LOO
  comparison against baselines.

## References

- Vehtari, Gelman, Gabry (2017). Practical Bayesian model evaluation using
  leave-one-out cross-validation and WAIC. *Statistics and Computing*.
- Bürkner, Gabry, Vehtari (2020). Approximate leave-future-out cross-validation
  for Bayesian time series models. *J. Stat. Comput. Simul.*
- Gelman (2006). Prior distributions for variance parameters in hierarchical
  models. *Bayesian Analysis*.

## Citation

```
@software{bayesianOU,
  author = {Gómez Julián, José Mauricio},
  title  = {bayesianOU: Bayesian Nonlinear Ornstein-Uhlenbeck Models},
  year   = {2024},
  note   = {ORCID: 0009-0000-2412-3150},
  url    = {https://github.com/IsadoreNabi/bayesianOU}
}
```

## License

MIT License
