# Porting notes

## Numerical translation

The optimizer was translated from the supplied package's f2c-derived C source
and checked against the corresponding Scilab fixed-form Fortran source. The
following numerical components are retained:

- packed LDL-transpose Hessian storage
- full-memory BFGS rank updates through `majour`
- the original search-direction triangular solves
- the original cubic interpolation and extrapolation line search
- automatic Hessian scaling from the `scale` vector
- mode-2 factorization of a supplied dense Hessian
- mode-3 reuse of an already factorized Hessian

The implementation is free-form Fortran 2018 with explicit typing, modules,
allocatable arrays, derived result/control types, IEEE finite checks, and
explicit callback interfaces.

## Hessian and restart distinction

The R wrapper returns `c.hess` as the packed lower triangle of the reconstructed
dense Hessian. The underlying mode-3 routine, however, expects the internal
packed LDL-transpose factor and skips factorization. Passing returned `c.hess`
with `restart=TRUE` therefore mixes two representations.

The Fortran API removes that ambiguity:

- `initial_hessian` accepts a dense Hessian and uses mode 2.
- `initial_factor` accepts `result%factor` and uses mode 3.
- `result%c_hess` is retained for R-output compatibility, not exact mode-3
  restart.

## Non-finite values

The translated API rejects an evaluation when either the objective or any
gradient component is non-finite. This is stricter and safer than the unusual
combined condition in the f2c-derived wrapper.

## Status reporting

The original R result does not expose a termination code. The Fortran result
adds statuses for convergence, limits, a non-descent direction, Hessian rank
loss, user cancellation, invalid input, and non-finite callback results.

## Omitted infrastructure

The following are R interface features rather than numerical algorithms and are
not translated:

- Rcpp `SEXP` evaluation
- R environments and assignment
- external function pointers registered with R
- dynamic-library registration
- R console formatting and S3-style list output
