# TMB — modern Fortran computational translation

This directory is a modern free-form Fortran/FPM translation of portable computational portions of
TMB 1.9.25 (Template Model Builder). It is intended to live as the top-level `TMB` directory in
`Fortran-from-R-packages`.

## Build

```text
fpm build
fpm test
fpm run --example simple
```

The maintained Fortran source is self-contained and uses no system BLAS/LAPACK link flags. It uses
one real kind, `dp`, defined from `real64` in `tmb_kinds` and re-exported by the public `tmb` module.

## Implemented computational API

The public `tmb` module currently exposes:

- scalar probability helpers: normal, exponential, Weibull, binomial, beta, gamma, lognormal, logistic,
  F, Student-t, skew-normal, multinomial, and sinh-asinh functions where implemented;
- numerically stable logit-binomial mass (`dbinom_robust`);
- matrix Kronecker product, sorting/order, and interval lookup;
- TMB-style smooth 2-D kernel interpolation with NaN data omission;
- Cholesky factorization and zero-mean multivariate-normal negative log density;
- standard-normal, stationary unit-variance AR(1), multivariate AR(1), and unstructured-correlation
  Gaussian negative log densities;
- matrix exponential by scaling/squaring with a convergent Taylor kernel;
- Romberg integration;
- deterministic central finite-difference gradient and Hessian helpers.

See `API_COVERAGE.md` for the explicit parity boundary.

## Important scope difference from upstream TMB

Upstream TMB is primarily a C++ automatic-differentiation and model-template framework built around
CppAD/TMBad, Eigen, R, and CHOLMOD. Reproducing that compiler/runtime architecture is not equivalent
to translating ordinary numerical routines. This package therefore translates reusable numerical
algorithms but does **not** claim drop-in parity with `MakeADFun`, C++ model templates, sparse AD
tapes, exact AD derivatives, CHOLMOD sparse factorization, or the R dynamic-library workflow.
Finite-difference derivative helpers are provided as a deterministic numerical facility, not as a
claim of exact-AD parity.

## Tests

`test/test_tmb.f90` checks distributions, Kronecker products, ordering, interpolation, Gaussian and
AR(1) negative log densities, matrix exponential, Romberg integration, and finite-difference
first/second derivatives with deterministic reference values.

## License and provenance

TMB is GPL-2. See `LICENSE` and `NOTICE`. Third-party code bundled with the original R package is
not vendored here.
