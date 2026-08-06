# tvgarch-modern-fortran

Modern Fortran computational translation of R package `tvgarch` 2.4.3.
The original package is by Susana Campos-Martins and Genaro Sucarrat and
is licensed GPL-2.0-or-later.

The library implements multiplicative time-varying GARCH models

```text
y(t) = sqrt(g(t) h(t)) z(t)
```

where `g(t)` is a smooth-transition long-term variance component and `h(t)`
is a GARCH/GJR-GARCH-X short-term component. Multivariate fitting estimates
margins equation by equation and supports constant or dynamic conditional
correlations.

## Implemented

- Logistic transition functions with one or several location parameters
- Multiple additive transition functions and all three upstream speed modes
- GARCH(p,q), GJR asymmetry, and nonnegative variance regressors
- Simulation with supplied or Gaussian innovations
- Maximization-by-parts estimation and optional joint refinement
- Robust sandwich covariance estimates for TV and GARCH blocks
- Fitted TV, GARCH, and total variance components
- Monte Carlo variance forecasts and empirical residual quantile paths
- Constant and Engle-style DCC correlations
- Multivariate simulation, marginal fitting, DCC fitting, and forecasting
- Lagged cross-variance spillover regressors
- Nonrobust and robust transition-order LM/TR2 tests
- Binary combination enumeration used by TV positivity constraints

R plotting, `zoo` indexes, formulas, S3 methods, printing, and LaTeX formatting
are intentionally omitted. Numeric arrays and typed result objects replace R
lists and classes.

## Build

GNU Fortran plus BLAS/LAPACK are required.

```text
make check
make optimized-check
```

An FPM manifest is included:

```text
fpm test
fpm run demo_tvgarch
```

FPM was unavailable in the validation environment; the manifest was parsed as
TOML and the equivalent source graph was tested with the Makefile.

## Main modules

- `tvgarch_transition`: TV transition functions and parameter handling
- `tvgarch_model`: univariate simulation, fitting, forecasting, and inference
- `tvgarch_multivariate`: CCC/DCC and multivariate/spillover workflows
- `tvgarch_tests`: transition-order specification tests
- `tvgarch`: umbrella public module

## Numerical notes

The supplied `garchx` translation is vendored and reused. Its bounded
Nelder-Mead optimizer replaces R's `constrOptim` and `nlminb`. Nonlinear TV
constraints are enforced by objective penalties and box bounds. Random streams
therefore differ from R, and optimizer endpoints are not claimed to be
bit-for-bit identical.

The GNU linker reports that optimizer procedure trampolines request an
executable stack. This originates from standard-conforming internal procedure
callbacks in the supplied `garchx` optimizer interface and the same callback
pattern in the TV/DCC layers. It does not affect the checked or optimized test
results, but consumers with a strict non-executable-stack policy may wish to
refactor the callback API to carry explicit context objects.
