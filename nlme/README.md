# nlme-fortran

A self-contained modern Fortran/FPM implementation of the numerical core of
R's `nlme` package for Gaussian linear and nonlinear mixed-effects models.

## Implemented computational areas

- Generalized least squares with ML or REML covariance estimation.
- Linear mixed-effects models with a general fixed-effects matrix `X`, random-
  effects matrix `Z`, and integer grouping vector.
- Nonlinear generalized least squares using typed model callbacks and
  Levenberg-Marquardt iterations.
- Nonlinear mixed-effects fitting by first-order iterative linearization.
- BLUP random effects, marginal and conditional fitted values and residuals,
  fixed-effect covariance, log likelihood, AIC, and BIC.
- AR(1), continuous AR(1), ARMA, compound-symmetry, exponential, Gaussian,
  linear, rational-quadratic, spherical, and unstructured correlations.
- Constant, fixed, group-specific, power, exponential, constant-plus-power,
  and constant-plus-proportional variance functions.
- Identity, diagonal, log-Cholesky, and compound-symmetry positive-definite
  random-effect covariance parameterizations.
- Residual ACF, empirical variograms, pooled standard deviations, grouped
  summaries, grouped linear/nonlinear fits, and mixed-model simulation.

## Design

R formulas, S3 classes, grouped-data objects, plotting, and lattice methods are
replaced by explicit arrays and Fortran derived types. The primary entry points
are:

```fortran
use nlme

call fit_gls(y, x, gls_fit, ...)
call fit_lme(y, x, z, group, lme_fit, ...)
call fit_gnls(model, y, xdata, theta0, gnls_fit, ...)
call fit_nlme(model, y, xdata, group, theta0, random_index, nlme_fit, ...)
```

All floating-point calculations use

```fortran
integer, parameter :: dp = kind(1.0d0)
```

The library has no external FPM dependencies.

## Build

```text
fpm build
fpm test
fpm run demo_nlme
```

Direct GNU Fortran validation scripts are included for systems without FPM.

## Important scope differences

This is a computational translation, not an implementation of R's formula/S3
runtime. Nested multi-level random effects are represented by constructing the
appropriate `Z` matrix and grouping vector explicitly. `fit_gnls` currently
uses supplied residual-covariance parameters as fixed; covariance optimization
is available in `fit_gls` and `fit_lme`. `fit_nlme` uses first-order
linearization rather than reproducing every detail of the original PNLS/native
optimizer stack. See `PORTING.md` and `TRANSLATION_COVERAGE.md`.
