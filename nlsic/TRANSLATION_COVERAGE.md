# Translation coverage

## Translated computational code

- `nlsic()` nonlinear constrained least squares:
  - projection of infeasible starts;
  - linear equalities and inequalities;
  - sequential constrained linearized least-squares steps;
  - backtracking globalization;
  - adaptive backtracking;
  - Jacobian reuse controls;
  - maximum step;
  - monotone-stop option;
  - optional least-norm step criterion;
  - numerical Jacobian fallback;
  - covariance, residual standard error, and t-based confidence widths;
  - parameter/direction/step/residual histories.
- `lsi()`.
- Current `lsi_ln()` least-norm constrained solver.
- `lsie_ln()`.
- `ldp()` using the Lawson-Hanson NNLS reduction.
- `ls_ln()` including rank-deficient and custom-norm cases.
- Matrix-right-hand-side `ls_ln()` behavior as `ls_ln_multi()`.
- `ls_ln_svd()` numerical role.
- `tls()`.
- `lsi_reg()`.
- `Nulla()` and `pnull()`.
- Numerical lower/upper bound conversion corresponding to `uplo2uco()`.

## R-only functionality omitted

`join()` and `g()` are string-formatting helpers and are not computational
numerical code.

`equa2vecmat()` parses and symbolically differentiates R expressions.  Fortran
callers supply the resulting numerical equality/inequality matrices directly.
The numerical lower/upper-bound part of `uplo2uco()` is provided by
`uplo_to_uco()`.

`lsi_lim()` is only an optional wrapper around the external R package
`limSolve` and is not duplicated.

The legacy, non-exported `lsi_ln_old()` and `ls_ln_old()` implementations are
retained in `original/` but are not duplicated beside the current algorithms.

## Numerical modernization

R's LAPACK QR/SVD operations are replaced by self-contained dense Gaussian
elimination and a symmetric Jacobi eigensolver used to construct
pseudoinverses/null spaces.  Rank cutoffs include a machine-precision floor to
avoid treating roundoff eigenvalues of `A^T A` as real singular directions.

Constrained rank-deficient least squares uses an active-set quadratic solve to
identify the active face, followed by an exact null-space least-norm
refinement.  This reproduces the mathematical role of the R implementation
without its QR-object attributes.

The R result attribute `aqr` is an optimization cache and is not exposed.
Internally the nonlinear solver still implements Jacobian reuse.

R `NA` values in nonlinear residuals correspond to non-finite residual entries
and are omitted from the nonlinear least-squares norm.  R-language exceptions
are represented by callback `ierr` values.

## Interface differences

The R `flsi` argument can inject an arbitrary R linear-solver function.  The
Fortran release exposes the translated `lsi`, `lsi_ln`, and `lsie_ln` routines
directly and selects least-norm nonlinear steps with
`control%least_norm_step`; an arbitrary procedure-valued `flsi` hook is not
yet part of the public `nlsic_solve` signature.

The R `lsi()` helper drops rows whose right-hand side is `NA`.  The standalone
Fortran linear routines expect finite numeric right-hand sides; non-finite
residual omission is implemented in the nonlinear `nlsic_solve` path where it
is algorithmically relevant.
