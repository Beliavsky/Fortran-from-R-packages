# Validation

## Environment

- Compiler: GNU Fortran (Debian 14.2.0-19) 14.2.0
- Language mode: Fortran 2018
- Linear algebra: system LAPACK and BLAS
- fpm: not installed
- Validation date: 2026-07-23

## Commands

Runtime-checked build and test workflow:

```sh
make check
```

Compiler flags:

```text
-std=f2018 -O0 -g -Wall -Wextra -Wimplicit-interface -Werror
-fcheck=all -fbacktrace
```

Optimized build and test workflow:

```sh
make release-check
```

Compiler flags:

```text
-std=f2018 -O2 -Wall -Wextra -Wimplicit-interface -Werror -fbacktrace
```

Both workflows link LAPACK and BLAS. Compiler warnings are treated as errors.

## Passing suites

```text
GPL-2.0-or-later source license checks passed.
Distribution tests passed.
Integration, grid, density, histogram, and binning tests passed.
Multivariate Normal, skew-Normal, skew-t, and skew-Cauchy fitting tests passed.
Original-name distribution and fitting compatibility tests passed.
```

The demonstration, integration example, and CSV fitter also ran successfully in
both workflows. The CSV fitter was exercised with Normal, skew-Normal,
fixed-degree skew-t, and skew-Cauchy modes.

## Numerical coverage

### Distribution tests

- Exact bivariate Normal origin identity
- Independent bivariate Normal and Student-t probabilities
- Cauchy equivalence to Student-t with one degree of freedom
- Normal and Student-t elliptical-density equivalence
- Positive logistic, Laplace, Kotz, and exponential-power densities
- Multivariate Normal and Student-t density positivity
- Multivariate Normal RNG means and covariance
- Multivariate Student-t RNG means and theoretical covariance
- Monte Carlo rectangular Normal probability against the quadrature-based
  bivariate Normal CDF
- Independent bivariate Normal equicoordinate quantile
- Zero-shape skew-Normal/Normal and skew-t/Student-t density equivalence
- Skew-Normal and skew-t simulation direction
- Skewed rectangular probability range and returned simulation error

### Fitting tests

- Multivariate Normal fitting
- One-dimensional skew-Normal fitting
- One-dimensional fixed- and free-degree skew-t fitting
- One-dimensional skew-Cauchy fitting
- Numerical Hessian and covariance allocation
- Positive-definite fitted scale matrices
- Two-dimensional skew-Normal fitting
- Two-dimensional fixed-degree skew-t fitting
- Two-dimensional skew-Cauchy fitting
- High-level method dispatch

### Integration and data-utility tests

- Original unit-square nine-point integration rule against the exact integral
  of `x*y`
- Adaptive 2-D integration against the exact integral of `exp(x+y)`
- Halton integration in two and three dimensions
- Cartesian-grid ordering and grid-data preservation
- KDE integral close to one
- Histogram count conservation
- Square- and hex-binning count conservation and positive nonempty counts

### Compatibility tests

Every original-name distribution wrapper and fit wrapper documented in
`API_MAP.md` is invoked. Compatibility aliases for `adapt`, `integrate2d`,
`delliptical2d`, `gridData`, `squareBinning`, and `hexBinning` are also
exercised.

## Accuracy and equivalence limits

- Bivariate Normal and Student-t CDFs use deterministic adaptive quadrature and
  are checked against analytical identities.
- Higher-dimensional and skewed rectangular probabilities are Monte Carlo
  estimates. Tests use returned standard errors and fixed seeds.
- The fitting tests establish executable likelihood, optimization, Hessian,
  covariance, fixed/free degrees-of-freedom, and multivariate paths. They do
  not claim equality to R's `sn` optimizer endpoints.
- R was not installed in the validation environment. Exact R random streams,
  optimizer traces, and external-package probability outputs were therefore
  not compared and are not claimed.
- The `fpm.toml` manifest is included but was not tested because fpm was not
  installed.
