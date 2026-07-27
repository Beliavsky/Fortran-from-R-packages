# garchx-modern-fortran

A modern Fortran computational translation of version 1.7 of the R package
`garchx` by Genaro Sucarrat.

The project translates numerical functionality. It does not reproduce R S3
classes, `zoo` indexes, plotting, printing, or LaTeX formatting.

## Implemented model

For selected positive lag sets, the conditional variance is

```text
sigma2(t) = omega
          + sum_j alpha(j) * y(t-arch_lag(j))^2
          + sum_j beta(j)  * sigma2(t-garch_lag(j))
          + sum_j gamma(j) * I[y(t-asym_lag(j)) < 0]
                           * y(t-asym_lag(j))^2
          + x(t)' delta
```

Lag sets may be sparse. For example, ARCH lags 1 and 4 can be included while
lags 2 and 3 are omitted.

## Implemented and tested

- Sparse-lag GARCH-X and asymmetric GARCH-X simulation
- Supplied or Gaussian standardized innovations
- Explicit/default backcast histories for simulation
- GARCH-X filtering and the package's Gaussian QML objective
- The zero-return objective mode (`objective.fun = 0` analogue)
- Bounded QML fitting using Nelder-Mead
- Exact recursive derivatives of conditional variances
- Numerical objective Hessians
- Ordinary Francq-Thieu covariance
- Robust covariance from the translated package formula
- HAC covariance with automatic or supplied bandwidth and Bartlett or supplied weights
- Gaussian log likelihood, fitted variance, and standardized residual output
- Bootstrap multi-step conditional-variance forecasts
- Empirical residual conditional-quantile paths
- Fixed-parameter and re-estimated refits
- Simulation-based ordinary asymptotic coefficient covariance
- Student-t confidence intervals
- Boundary-null one-sided t tests
- Simulation-based boundary Wald tests
- Multivariate Normal simulation
- Vector and matrix lag/difference helpers
- Dated CSV fitting application

See `API_MAP.md` for the mapping from R functions to Fortran procedures and
`VALIDATION.md` for the executed checks.

## Build and test

GNU Fortran, LAPACK, and BLAS are required.

```text
make check
make optimized-check
```

The runtime-checked build uses bounds, allocation, and other runtime checks and
treats warnings as errors. `fpm.toml` is also included:

```text
fpm test
fpm run demo_garchx
fpm run fit_csv -- data/example.csv 1 1 1 hac
```

`fpm` was not available in the validation environment, so only the Makefile
commands are claimed as tested.

## CSV application

The CSV format is:

```text
Date,y[,xreg1,xreg2,...]
```

Examples:

```text
build/debug/fit_csv data/example.csv 1 1 - ordinary
build/debug/fit_csv data/example.csv 1,4 1,3 2 hac
```

A dash means no selected lags. Any columns after `y` are treated as covariates.
The parameter order is intercept, selected ARCH coefficients, selected GARCH
coefficients, selected asymmetry coefficients, and covariate coefficients.

## Important differences from R

- Bounded Nelder-Mead replaces R's `nlminb`; exact optimizer endpoints are not claimed.
- The Fortran random stream differs from R's random stream.
- Numeric arrays replace formulas, `zoo` objects, and name-based extraction.
- Bounds must currently be supplied as full vectors in the library API.
- Missing values must be removed before calling the numerical API.
- `garchxAvar` in the original package implements ordinary covariance only;
  this translation likewise does not claim robust or HAC `garchxAvar` modes.
- Exact output equivalence with R is not claimed because R was unavailable in
  the validation environment.

## License

The original package declares `GPL (>= 2)`. This translation is distributed as
`GPL-2.0-or-later`. `LICENSE` contains GNU GPL version 2, and each Fortran source
contains an SPDX identifier, the version-2-or-later notice, and original-package
attribution.
