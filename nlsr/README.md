# nlsr-fortran

Modern Fortran translation of the computational core of the R package **nlsr 2026.4.29** by John C. Nash, Duncan Murdoch, Fernando Miguez, and Arkajyoti Bhattacharjee.

The central routine is `nlfb`, Nash's stabilized nonlinear least-squares solver.  It works directly with residual and optional residual-Jacobian callbacks and supports parameter bounds, exact masks (`lower == upper`), observation weights, callback-generated weights, Marquardt stabilization, backtracking, and the residual-orthogonality (`roff`) convergence test used by upstream nlsr.

## Build

```text
fpm build
fpm test
```

No BLAS, LAPACK, R, or external Fortran library is required.

## Minimal use

```fortran
use nlsr

type(nlsr_result) :: fit
real(dp) :: start(2)
start = [1.0_dp, 1.0_dp]
call nlfb(start, nres, residual, fit, jacfn=jacobian)
```

Callbacks have explicit interfaces:

```fortran
subroutine residual(par, r, ierr)
    use nlsr, only : dp
    real(dp), intent(in) :: par(:)
    real(dp), intent(out) :: r(:)
    integer, intent(out) :: ierr
end subroutine

subroutine jacobian(par, jac, ierr)
    use nlsr, only : dp
    real(dp), intent(in) :: par(:)
    real(dp), intent(out) :: jac(:,:)
    integer, intent(out) :: ierr
end subroutine
```

When `jacfn` is omitted, the `nlsr_control%jacobian_method` setting selects forward, backward, central, or Richardson finite differences.

See `API.md`, `TRANSLATION_COVERAGE.md`, and the `example/` directory.
