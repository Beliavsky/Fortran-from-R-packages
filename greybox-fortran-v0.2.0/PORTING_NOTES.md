# Porting notes

## General design

The upstream package mixes numerical modeling with R formula/model-frame/S3/plotting infrastructure. The Fortran port exposes arrays, matrices, scalars and derived result types. It intentionally does not emulate formulas, model frames, `zoo` indices, graphics, `xtable` or `texreg` output.

## Optimisation

Upstream `alm` can use `nloptr`, `optim`, and Hessians from `pracma`. The Fortran port uses a self-contained deterministic coordinate pattern search plus finite-difference Hessians; ordinary MSE regression uses direct least squares. A `pracma` translation was available during development, but only a small Hessian routine was required and is implemented locally.

This means difficult non-convex likelihoods can settle at a different local optimum from R even when the likelihood and parameterization are the same.

## ALM response and latent scales

`alm_model%fitted` and `alm_predict` are response-scale values. Likelihood evaluation for lognormal/log-Laplace/log-S/log-generalized-normal, Box-Cox-normal and logit-normal families uses the corresponding latent parameterization.

## Beta regression

Upstream `alm(..., distribution="dbeta")` uses two coefficient blocks:

- `shape1 = exp(X beta1)`
- `shape2 = exp(X beta2)`
- fitted mean `shape1/(shape1+shape2)`

The Fortran port follows that parameterization. `alm_model%beta` contains `beta1`, `scale_beta` contains `beta2`, and the Hessian covariance is exposed in three blocks: `vcov`, `vcov_scale`, and `vcov_cross`. Boundary observations at 0/1 receive the same small inward correction idea used upstream.

`calm_fit` also averages both coefficient blocks for beta models.

## Occurrence / hurdle modeling

The R function folds occurrence modeling into a single large `alm` object. Fortran exposes the same numerical composition through `alm_fit_occurrence`: a Bernoulli logistic/probit occurrence model is estimated on all observations, the positive-demand distribution is estimated on positive observations, and the combined point likelihood and response-scale fitted mean are returned in `alm_occurrence_model`.

This separate API is intentional; it avoids reproducing R's formula/object dispatch while retaining the statistical model.

## ARIMA-error orders

`alm_fit_arima_errors` implements the numerical role of `orders=c(p,d,q)` using conditional iteration:

1. difference `y` and the design matrix `d` times;
2. augment the design with `p` lagged transformed responses and `q` lagged residuals;
3. refit the selected ALM distribution;
4. update residuals and repeat.

Upstream optimizes the recursive likelihood jointly inside `alm`; therefore AR/MA estimates can differ slightly. The Fortran method is a stable conditional-likelihood implementation and preserves the same model structure without requiring `nloptr`.

## Robust and penalized losses

The v0.2 implementation includes:

- `MSE`, `MAE`, `HAM`
- `LASSO`: `(1-lambda) MSE + lambda sum(abs(beta))`
- `RIDGE`: `(1-lambda) MSE + lambda ||beta_without_intercept||_2`
- `ROLE`: trimmed mean point log likelihood
- `QUALE`: selected sample quantile of point log likelihood

ROLE uses the same two-sided trimming concept as upstream and QUALE uses R type-7-compatible linear quantile interpolation. For LASSO/RIDGE, point likelihoods are still available but overall likelihood-based IC values are not treated as the fitting objective, matching the conceptual distinction upstream makes for penalized-loss fits.

## `lmDynamic`

`lm_dynamic_fit` enumerates the same model subsets as exhaustive `lmDynamic`, obtains point IC values for every observation/model pair, forms row-wise Akaike/Bayesian weights, logit-transforms the weights, smooths them with a locally weighted linear tricube smoother with robust iterations, inverse-logit transforms and renormalizes them, and returns dynamic coefficients and variable importance.

When no span is supplied, the Fortran routine chooses among spans 0.2 through 1.0 using the same weighted point-IC objective idea. Upstream continuously optimizes the LOWESS span through `nloptr`; this grid search is deterministic and avoids the external optimizer.

## `rmcb`

The Tukey branch uses row ranks, the Friedman statistic, and an infinite-degrees-of-freedom studentized-range quantile. The latter is evaluated numerically from the exact range distribution of independent standard normals. The `dnorm` branch uses the rank regression/F test formulation. Other distributions use `alm_fit` on the original values and a likelihood-ratio comparison against the intercept-only model.

## `dsrboot`

The port covers additive/multiplicative, parametric/nonparametric, scaling and intermittent-demand branches. It preserves the upstream sorted-series perturbation/reordering/centering construction. Upstream uses `supsmu` for several smoothing steps; the Fortran version uses a deterministic local moving smoother in `dsrboot` because R's `supsmu` implementation is not part of the greybox source tree. The bootstrap law and zero/intermittent reconstruction remain the same workflow, but individual simulated paths will not match R seed-for-seed.

## `aid` / `aidCat`

The Fortran demand identifier fits the same principal candidate families (regular/intermittent, fractional/count, smooth/lumpy), uses occurrence modeling for intermittent candidates, returns the six upstream category names, and records new/obsolete/stockout anomalies. `aid_cat` aggregates the same 2-by-3 category table and anomaly counts.

The source R implementation relies heavily on `supsmu`/`lowess` and neighbor-based stockout heuristics. The Fortran routine uses its native smoother and a fitted-geometric tail threshold for internal zero runs, so these classifiers are numerically equivalent in intent but not expected to make identical calls on borderline series.

## Distribution details retained from v0.1

- generalized-normal large-shape limiting-uniform handling;
- rectified-normal atom at zero;
- Box-Cox transformation/Jacobian parameterization;
- response-scale predictions for transformed families.

## Information criteria

AIC, AICc, BIC and BICc follow `greybox`. Point criteria preserve the upstream convention in which each observation receives the full parameter penalty, so the mean of `point_aic(model)` equals model AIC and similarly for BIC.

## Remaining worthwhile numerical gaps

After v0.2, the remaining differences are narrower:

1. joint `nloptr` optimization of recursive ARIMA-error ALM rather than conditional iteration;
2. exact `supsmu` behavior inside `aid` and `dsrboot`;
3. custom user-supplied ALM loss callbacks;
4. some scale/shape submodel combinations and rarely used distribution-specific Fisher-information shortcuts;
5. R's parallel execution paths.

These are lower priority than the numerical parity targets closed in v0.2. Plotting and R presentation/object infrastructure remain intentionally out of scope.

## Licensing

Current upstream `DESCRIPTION` declares `LGPL-2.1`. An older package-level documentation comment says GPL-2. The current package metadata is treated as authoritative; historical source text remains under `orig/` for provenance.
