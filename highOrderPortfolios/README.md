# highOrderPortfolios-fortran

Modern Fortran translation of the computational core of the R package
`highOrderPortfolios` 0.1.1. The project is an FPM library and has no external
Fortran dependencies.

The library designs long-only portfolios using mean, variance, co-skewness,
and co-kurtosis. It supports sample-moment and generalized hyperbolic skew-t
parameterizations, as well as portfolio tilting under a tracking-error bound.

## Implemented public routines

- `estimate_sample_moments`
- `estimate_skew_t`
- `eval_portfolio_moments`
- `design_mvsk_portfolio_via_sample_moments`
- `design_mvsk_portfolio_via_skew_t`
- `design_mvsktilting_portfolio_via_sample_moments`

R names containing capitals are written in lower-case Fortran style. The
Fortran compiler treats names case-insensitively, so the recognizable spelling
may still be used in source code.

## Design choices

`sample_moments` stores the centered observations and evaluates third- and
fourth-order portfolio moments directly. This requires O(Tp + p^2) storage,
rather than materializing an O(p^4) cokurtosis array. Explicit co-skewness and
cokurtosis tensors can be requested with `store_tensors=.true.`.

The skew-t estimator is adapted from the previously translated
`fitHeavyTail-fortran` project and is included in this source tree. The R
optimization dependencies (`quadprog`, `ECOSolveR`, `lpSolveAPI`, and `nloptr`)
are replaced by self-contained simplex projection, convex quadratic projected
solvers, successive convex approximation, accelerated projected-gradient, and
tracking-feasible tilting routines.

## Build

```text
fpm build
fpm test
fpm run
```

Without FPM, run `scripts/test_gfortran.sh` on Unix-like systems or
`scripts/test_gfortran.bat` on Windows.

## Licensing

This translation is licensed under GPL-3.0-only, matching the upstream
package and the reused `fitHeavyTail-fortran` numerical code. See `LICENSE` and
`NOTICE.md`.
