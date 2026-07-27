# lgarch-modern-fortran

A modern Fortran translation of the computational parts of the CRAN package
`lgarch` 0.7 by Genaro Sucarrat.

The translation preserves the original package's GPL-2 license as
`GPL-2.0-only`. Every Fortran source file contains an SPDX identifier and an
explicit GPL version 2-only notice.

## Implemented numerical functionality

### General utilities

- Lagging vectors and matrices with optional padding.
- Differencing vectors and matrices with optional padding.
- Multivariate Normal random generation using a Cholesky factor.
- LAPACK/BLAS-backed solves, inverses, Cholesky factors, log determinants,
  eigenvalue stability checks, covariance matrices, and correlation matrices.

### Univariate log-GARCH

- `lgarchSim` equivalent for arbitrary simulation orders `p` and `q`.
- User-supplied or Normal innovations.
- Optional contemporaneous forcing series, initial log-variance values, and
  initial log-squared-innovation values.
- AR-polynomial stability check.
- ARMA recursion with the original zero-return treatment: a zero observation
  is replaced by its conditional fitted log-square and contributes a zero ARMA
  residual.
- Least squares, Gaussian QML, and centred exponential chi-squared QML.
- Optional mean correction.
- Exogenous regressors in estimation.
- Custom starting values, lower bounds, upper bounds, optimization tolerance,
  iteration limit, and invalid-objective penalty.
- Numerical Hessian and ARMA covariance matrix.
- Transformation from ARMA coefficients to log-GARCH coefficients.
- Fitted conditional standard deviations and log variances.
- Standardized log-GARCH residuals and ARMA residuals.
- Log-GARCH Gaussian log likelihood, ARMA log likelihood, and ARMA RSS.

The original `lgarch()` estimator restricts `arch` and `garch` to orders zero
or one, with `garch <= arch`. The Fortran estimator intentionally preserves
that restriction. The simulation routine supports higher orders, as does the
original `lgarchSim()` routine.

### Multivariate CCC log-GARCH

- Multivariate Normal innovation generation with a supplied covariance matrix.
- CCC multivariate log-GARCH(1,1) simulation.
- Optional supplied innovations, equation-specific forcing series, and
  backcast values.
- Stability checking from the eigenvalues of `ARCH + GARCH`.
- VARMA recursion with elementwise zero-return imputation.
- Gaussian VARMA log likelihood with a full residual covariance matrix.
- Bounded multivariate estimation for order-zero and order-one specifications.
- Custom starting values, bounds, optimizer controls, and objective penalty.
- Parameter transformation from VARMA to multivariate log-GARCH form.
- Fitted conditional standard deviations, log variances, standardized
  residuals, and VARMA residuals.
- VARMA Hessian and covariance matrix.
- Gaussian CCC log-GARCH likelihood based on the fitted residual correlation.

## Result types

`lgarch_fit_result` contains the computational information otherwise exposed
through the R `coef`, `fitted`, `logLik`, `residuals`, `rss`, and `vcov`
methods. `mlgarch_fit_result` provides the corresponding multivariate values.
No R-style runtime class system is reproduced.

For the transformed log-GARCH covariance matrices, the dynamic and exogenous
coefficient blocks are obtained by a delta-method transformation holding the
estimated `Elnz2` values fixed. The separate `Elnz2` variances follow the
package formulas. Cross-covariances involving `Elnz2`, and unidentified
intercept terms, remain IEEE NaN rather than being presented as estimated.

## Building and testing

GNU Fortran, LAPACK, and BLAS are required.

```sh
make check
```

The checked build uses Fortran 2018, full runtime checking, backtraces, and
warnings as errors. An optimized validation can be run with:

```sh
make clean
make FFLAGS="-std=f2018 -O2 -Wall -Wextra -Wimplicit-interface -Werror -ffree-line-length-none -fbacktrace" check
```

A valid `fpm.toml` is included. The validation environment did not contain
`fpm`, so `fpm build` is not claimed as tested.

## Programs

```sh
build/demo_lgarch
build/xreg_example
build/fit_csv data/univariate.csv univariate ls 1 1
build/fit_csv data/univariate.csv univariate ml 1 1
build/fit_csv data/univariate.csv univariate cex2 1 1
build/fit_csv data/multivariate.csv multivariate 1 1
```

The CSV format is `Date,Series1,...`. The first field is read as a date label
and ignored numerically.

## Differences from R

- The Fortran optimizer is bounded Nelder-Mead rather than R's `nlminb`.
- Numerical Hessians use central finite differences rather than `optimHess`.
- Random-number streams are native Fortran streams and do not reproduce R's
  RNG sequence.
- R formulas, `zoo` indexes, S3 classes, printing methods, summaries, and
  plotting are excluded.
- Exact coefficient-by-coefficient equality with R optimization output is not
  claimed.
- The original source contains commented-out skew-GED and skew-t generators,
  asymmetry placeholders, and unimplemented prediction methods. They are not
  exported computational features of `lgarch` 0.7 and are not claimed here.

See `API_MAP.md`, `ORIGIN.md`, and `VALIDATION.md` for details.
