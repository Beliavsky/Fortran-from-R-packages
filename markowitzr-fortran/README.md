# markowitzr-fortran

A self-contained modern Fortran/FPM translation of the computational routines
in R package `MarkowitzR` 1.0.2.0002.

The library estimates unified second moments, applies the delta method to their
inverses, and estimates conditional or unconditional Markowitz coefficients
under optional subspace and hedging constraints.

## Build

```text
fpm build
fpm test
fpm run markowitzr_demo
fpm run --example basic_inference
fpm run --example conditional_constraints
```

A direct GNU Fortran validation script is also included:

```text
scripts/validate.sh
scripts\validate.bat
```

## Main interfaces

```fortran
use markowitzr, only: dp, theta_result, markowitz_result
use markowitzr, only: theta_vcov, itheta_vcov, mp_vcov
```

### Unified second moment

```fortran
result = theta_vcov(x)
```

For an `n x p` return matrix and the default fitted intercept, `result%mu`
contains

```text
vech(E([1,x]' [1,x]))
```

in lower-triangular, column-major order. `result%covariance` is the estimated
covariance matrix of that packed sample moment.

Available covariance estimators are:

```fortran
covariance_empirical
covariance_normal
covariance_hac
```

A user procedure matching `moment_covariance_callback` may be supplied instead.
The callback receives the observation-by-moment matrix and returns a covariance
matrix for its column means.

### Inverse unified second moment

```fortran
result = itheta_vcov(x)
```

This returns `vech(inv(Theta))` and its analytic delta-method covariance. With
an intercept, the lower-left block contains the negative unscaled Markowitz
portfolio and the lower-right block contains the precision matrix.

### Markowitz coefficient inference

```fortran
fit = mp_vcov(x, feat=features, jmat=jmat, gmat=gmat)
```

Important fields are:

```fortran
fit%w              ! p x ff Markowitz coefficient matrix
fit%w_covariance   ! covariance of vec(w), column-major
fit%mu             ! vech of the projected inverse moment
fit%covariance     ! covariance of fit%mu
fit%w_indices      ! positions of vec(w) within fit%mu
fit%n
fit%p
fit%ff
fit%status
fit%message
```

`jmat` constrains portfolios to its row space. `gmat` imposes zero-covariance
hedging constraints. When both are supplied, the rows of `gmat` must lie in the
row space of `jmat`.

## Weight compatibility

The upstream R implementation applies `weights` only to feature columns, even
though its documentation says returns and features are multiplied. The default
Fortran mode preserves that behavior:

```fortran
weight_mode=weights_upstream
```

The documented all-column transformation is available explicitly:

```fortran
weight_mode=weights_all_columns
```

This distinction is recorded in `PORTING_NOTES.md`.

## Numerical design

- Fortran 2018 modules and typed result structures
- `dp = kind(1.0d0)`
- no external linear algebra or statistics dependencies
- partial-pivoting matrix inversion
- analytic duplication/Kronecker delta-method Jacobians
- exact covariance symmetrization
- complete-case row handling
- native empirical, Gaussian, and Bartlett-HAC covariance estimators
- procedure callback support for custom covariance estimators

## Scope

All three functions exported by the R package are translated:

- `theta_vcov`
- `itheta_vcov`
- `mp_vcov`

The small internal matrix helpers required by those routines are also included.
R formula objects, `lm` objects, package documentation machinery, and optional
`sandwich` adapters are not compiled. Their numerical role is replaced by
native estimators and the callback interface.

## License

`LGPL-3.0-or-later`, matching the original source headers. See `LICENSE`,
`COPYING.GPL`, and `NOTICE`.
