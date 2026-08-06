# Porting notes

## Model parameterization

For both models, the conditional mean is

```text
mu_i = exp(a_i^T beta).
```

The exponential negative log-likelihood follows the upstream C++ scaling and
is the unnormalized sum of `eta_i + y_i exp(-eta_i)`. The fixed-shape Gamma
negative log-likelihood follows the R implementation and is averaged over
observations.

## Gamma shape estimation

The R code passes the shape directly to unconstrained BFGS. The Fortran port
optimizes `log(shape)`, which guarantees positivity. It uses an analytic
shape gradient, a self-contained digamma approximation, inverse-BFGS updates,
and Armijo backtracking. The estimated shape is then held fixed during Gamma
elastic-net cross-validation, matching the package workflow.

## Elastic-net proximal map

The correct proximal update for a penalized coefficient is

```text
soft(v, step*lambda*alpha) / (1 + step*lambda*(1-alpha)).
```

The upstream C++ code instead divides only the threshold by the denominator.
The corrected map is the default. Set `source_proximal=.true.` to reproduce
the C++ expression.

The R Gamma implementation penalizes all coefficients. The C++ exponential
implementation leaves the first coefficient unpenalized when an intercept is
present. `penalize_intercept` makes this choice explicit. The high-level Gamma
wrapper adopts the R behavior unless the caller supplies options.

## Cross-validation corrections

Two C++ indexing expressions cannot produce valid folds:

1. `(k_fold-1)/k_fold*num_obs` uses integer division and is zero for normal
   `k_fold > 1`;
2. the test-index vector is written using original observation indices rather
   than positions starting at zero.

The port uses deterministic Fisher-Yates permutations and balanced fold IDs.
`seed` and `k_fold_iter` control repeatability.

The source chooses lambda by Euclidean prediction error on `exp(A beta)`.
`cv_metric='source'`, `'prediction'`, or `'rmse'` reproduces response-scale
selection. The default `cv_metric='nll'` uses held-out negative log-likelihood.

## Lambda grid

The Fortran grid is descending from `lambda_max` to
`lambda_max*min_lambda_ratio`, enabling warm starts. `min_lambda_absolute`
can request an absolute lower endpoint, including the upstream exponential
choice of `0.001`.

For pure ridge (`alpha=0`), no finite lambda makes all penalized coefficients
exactly zero. The implementation uses a small alpha floor only to establish a
practical grid scale.

## Omitted infrastructure

Rcpp registration, R input-type checks, formula/model objects, and vignette
runtime dependencies are not part of the numerical library. The trivial
`MyClass` addition and scalar-multiplication examples are also omitted.
