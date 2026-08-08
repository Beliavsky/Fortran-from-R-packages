# API

Main module: `marqlevalg`

## Types

### `mla_control`

- `maxiter = 500`
- `epsa = 1e-4`
- `epsb = 1e-4`
- `epsd = 1e-4`
- `minimize = .true.`
- `blinding = .true.`
- `multiple_try = 25`
- `partial_h(:)` optional 1-based parameter indices excluded from the RDM
  information matrix when the full matrix is singular.

### `mla_result`

- `par(:)` final parameters
- `fn_value` objective on the user's original scale
- `grad(:)` final gradient on the user's original objective scale
- `vcov(:,:)` inverse information matrix when invertible; otherwise the final
  information matrix
- `iterations`
- `ier` information-matrix inversion status
- `istop`: 1 full convergence, 2 iteration limit, 3 partial-Hessian
  convergence, 4 numerical failure
- `ca`, `cb`, `rdm`

## Optimizer

```fortran
call marqlev_optimize(x0, fn, result [, control])
call marqlev_optimize(x0, fn, gr, result [, control])
call marqlev_optimize(x0, fn, gr, hess, result [, control])
```

Callback interfaces are explicit and assumed-shape:

```fortran
function fn(x) result(f)
    real(dp), intent(in) :: x(:)
    real(dp) :: f
end function

subroutine gr(x, g)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: g(:)
end subroutine

subroutine hess(x, h)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: h(:, :)
end subroutine
```

## Derivatives

```fortran
call deriva(x, fn, f, info, grad)
call deriva_grad(x, gr, info)
```

These preserve the upstream finite-difference step rule
`max(1e-7, 1e-4*abs(x(i)))` and its forward mixed-second-derivative formula.

## Mixed model utilities

Module `mla_lmm`:

```fortran
ll = loglik_lmm(b, y, x, ni)
call grad_lmm(b, y, x, ni, grad)
```

They implement the upstream random-intercept Gaussian linear mixed model.
