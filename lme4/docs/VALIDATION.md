# Validation

## Automated tests

The runtime-checked test program covers:

- covariance and standard-deviation-correlation round trips;
- unstructured, diagonal, compound-symmetry, and AR(1) covariance construction;
- standard-normal Gauss-Hermite moments;
- random-intercept, random-slope, and multiple-term Gaussian LMMs;
- equivalence of dense and Woodbury/PLS LMM fixed effects, residual scale, and likelihood;
- binomial-logit, Poisson-log, Gamma-log, inverse-Gaussian-log, and negative-binomial-log GLMMs;
- negative-binomial size profiling;
- scalar adaptive quadrature for binomial and Poisson models;
- correlated two-dimensional adaptive quadrature for a binomial random-intercept/random-slope model;
- custom binomial-probit family callbacks;
- a Gaussian nonlinear mixed growth model;
- Wald and profile-likelihood fixed-effect intervals;
- parametric bootstrap and percentile intervals;
- group-deletion influence refits;
- likelihood-ratio probabilities;
- grouped `lmList`-style regressions and prediction;
- deterministic simulation for implemented response families.

## Dense versus Woodbury LMM

A deterministic random-intercept model is fitted by both `fit_lmm` and `fit_lmm_pls`. The automated tolerances require:

| Quantity | Required maximum difference |
|---|---:|
| Fixed effects | `1e-6` |
| Deviance | `1e-5` |
| Residual scale | `1e-6` |

The observed optimized-build differences are approximately `4e-11` for fixed effects and `1.3e-9` for deviance.

## Independent LMM comparison

A deterministic 60-observation random-intercept model was fitted by this port and Python `statsmodels` 0.14.6 using ML.

| Quantity | Maximum absolute difference |
|---|---:|
| Fixed effects | `2.9e-14` |
| Random-intercept variance | `2.7e-05` |
| Residual variance | `1.1e-07` |
| Log likelihood | `1.2e-08` |

## Extended-family checks

Synthetic data are generated with fixed seeds and refitted. Tests require convergence, valid fitted-value ranges, retention of fixed dispersion, and recovery of the simulated slope within Monte Carlo tolerances. The negative-binomial profile must improve on a boundary-size fit.

## Adaptive quadrature checks

Scalar AGHQ is checked for binomial and Poisson models. Multidimensional AGHQ is checked on a correlated random-intercept/random-slope binomial model for convergence, positive variances, valid response-scale fitted values, and correct covariance dimension.

A separate deterministic 15-node scalar AGHQ comparison against direct infinite-interval numerical integration gives:

| Quantity | Maximum absolute difference |
|---|---:|
| Fixed effects and random-effect SD | `1.6e-08` |
| Marginal log likelihood | `1.2e-13` |

The independent calculation is in `validation/aghq_validation.py`.

## Nonlinear model check

A deterministic exponential-growth model with grouped random intercepts is refitted from displaced starting values. The recovered fixed parameters are approximately `1.39546` and `0.55141` for generating values `1.4` and `0.55`.

## Inference checks

The profile-likelihood routine is exercised on a Gaussian random-intercept model and must bracket both ML fixed-effect estimates. Parametric bootstrap requires successful refits, influence diagnostics require successful deletion fits for nearly all groups, and the chi-square implementation reproduces the one-degree-of-freedom tail probability for a likelihood-ratio statistic of four (`0.0455003`).
