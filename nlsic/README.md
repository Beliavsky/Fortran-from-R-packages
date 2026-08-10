# nlsic-fortran

Modern Fortran/FPM translation of the computational code in the R package
`nlsic` 1.2.0 (Serguei Sokol).

The package solves nonlinear least-squares problems with linear equality and
inequality constraints.  The nonlinear iteration follows `nlsic()`:
sequential constrained linear least-squares steps, feasible-start projection,
and backtracking globalization.

## Main API

```fortran
use nlsic

type(nlsic_control) :: control
type(nlsic_result)  :: fit

call nlsic_solve(par0, nres, residual, fit, &
   jacobian=jac, u=u, co=co, e=e, eco=eco, control=control)
```

`residual` has the interface

```fortran
subroutine residual(par, r, ierr)
   real(dp), intent(in)  :: par(:)
   real(dp), intent(out) :: r(:)
   integer, intent(out)  :: ierr
end subroutine
```

and an analytical Jacobian callback, when supplied, has the interface

```fortran
subroutine jacobian(par, r, j, ierr)
   real(dp), intent(in)  :: par(:)
   real(dp), intent(out) :: r(:), j(:,:)
   integer, intent(out)  :: ierr
end subroutine
```

If no Jacobian callback is supplied, a central finite-difference Jacobian is
computed internally.

## Linear numerical routines

The translated package also exposes the numerical helpers that are exported by
the R package:

- `ldp` -- least-distance programming, `U*x >= co`.
- `lsi` -- full-rank least squares with inequalities.
- `lsi_ln` -- constrained least squares with least-norm tie breaking.
- `lsie_ln` -- least-norm least squares with equalities and inequalities.
- `ls_ln` -- least-norm least squares.
- `ls_ln_multi` -- multiple right-hand sides.
- `ls_ln_svd` -- SVD/pseudoinverse-style least-norm solution.
- `lsi_reg` -- rank-deficient regularized constrained least squares.
- `tls` -- total least squares.
- `nulla` -- null-space basis (vector or matrix input).
- `pnull` -- particular solution plus null-space basis.
- `uplo_to_uco` -- convert numeric lower/upper bounds to `U*x >= co`.

The implementation is standalone and requires no BLAS/LAPACK, `nnls`,
`numDeriv`, `limSolve`, `dotty`, or `glue` dependency.

## Build

```text
fpm build
fpm test
```

GNU Fortran strict validation scripts are in `scripts/`.

## License

The upstream package is GPL-2.0-only.  The original source archive is retained
under `original/nlsic-master/` and `COPYING` contains GPL version 2.
