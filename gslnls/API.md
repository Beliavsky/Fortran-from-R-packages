# Fortran API

All public symbols are re-exported by:

```fortran
use gslnls
```

## Main types

### `type(nls_control)`

Important fields:

- `maxiter`
- `algorithm`: `NLS_LM`, `NLS_LMACCEL`, `NLS_DOGLEG`, `NLS_DDOGLEG`,
  `NLS_SUBSPACE2D`, or `NLS_CGST`
- `scale`: `NLS_SCALE_MORE`, `NLS_SCALE_LEVENBERG`, or
  `NLS_SCALE_MARQUARDT`
- `solver`: `NLS_SOLVER_QR`, `NLS_SOLVER_CHOLESKY`, or `NLS_SOLVER_SVD`
- `fdtype`: `NLS_FD_FORWARD` or `NLS_FD_CENTER`
- `factor_up`, `factor_down`, `avmax`
- `h_df`, `h_fvv`
- `xtol`, `ftol`, `gtol`
- multi-start fields `mstart_n`, `mstart_p`, `mstart_q`, `mstart_r`,
  `mstart_s`, `mstart_tol`, `mstart_maxiter`, `mstart_maxstart`,
  `mstart_minsp`
- robust-fit fields `irls_maxiter`, `irls_xtol`
- `store_trace`

Defaults mirror the R package where practical.

### `type(nls_loss)`

Construct with:

```fortran
loss = default_loss(LOSS_HUBER)
```

Available constants are `LOSS_DEFAULT`, `LOSS_HUBER`, `LOSS_BARRON`,
`LOSS_BISQUARE`, `LOSS_WELSH`, `LOSS_OPTIMAL`, `LOSS_HAMPEL`, `LOSS_GGW`, and
`LOSS_LQQ`.

### `type(nls_result)`

Contains:

- `par`, `fitted`, `residual`, `weighted_residual`
- final `jacobian`, `covariance`
- robust `irls_weights`, `irls_psi`, `irls_dpsi`
- optional `par_trace`, `ssr_trace`
- `ssr`, `sigma`, `irls_sigma`, `gradient_inf`
- iteration/evaluation counters
- numerical `rank`
- `status`, `info`, `converged`

Status constants are `NLS_SUCCESS`, `NLS_MAXITER`, `NLS_NO_PROGRESS`,
`NLS_BAD_FUNCTION`, `NLS_BAD_INPUT`, and `NLS_SINGULAR`.

## Callback interfaces

```fortran
subroutine model(par, yhat, ierr)
  real(dp), intent(in) :: par(:)
  real(dp), intent(out) :: yhat(:)
  integer, intent(out) :: ierr
end subroutine
```

An analytic Jacobian has the form:

```fortran
subroutine jacobian(par, jac, ierr)
  real(dp), intent(in) :: par(:)
  real(dp), intent(out) :: jac(:,:)
  integer, intent(out) :: ierr
end subroutine
```

A geodesic directional second derivative has the form:

```fortran
subroutine fvv(par, v, value, ierr)
  real(dp), intent(in) :: par(:), v(:)
  real(dp), intent(out) :: value(:)
  integer, intent(out) :: ierr
end subroutine
```

The matrix-free large-system interface uses:

```fortran
subroutine jacobian_operator(par, transpose_j, u, v, ierr)
  real(dp), intent(in) :: par(:), u(:)
  logical, intent(in) :: transpose_j
  real(dp), intent(out) :: v(:)
  integer, intent(out) :: ierr
end subroutine
```

When `transpose_j=.false.`, return `J*u`; when true, return `transpose(J)*u`.

## Fitting routines

### `fit_nls`

```fortran
call fit_nls(model, y, start, result, control, jac, fvv, &
             lower, upper, weights, weight_matrix, loss)
```

All arguments after `result` are optional. `weights` and `weight_matrix` are
mutually exclusive.

### `fit_nls_large`

Dense large-system convenience entry point. It selects `NLS_CGST` by default
and otherwise shares the `fit_nls` representation.

### `fit_nls_large_operator`

```fortran
call fit_nls_large_operator(model, jop, y, start, result, control, weights)
```

This avoids forming a dense Jacobian and applies Steihaug CG through Jacobian
operator products.

### `fit_nls_multistart`

```fortran
real(dp) :: ranges(2,p)
type(multistart_result) :: ms
call fit_nls_multistart(model, y, ranges, ms, control, jac=jacobian)
```

`ranges(1,:)` and `ranges(2,:)` are lower and upper starting ranges.

## Numerical derivative helpers

- `numerical_jacobian`
- `numerical_fvv`

## Inference helpers

- `nls_hatvalues`
- `nls_cooks_distance`
- `nls_loglik`
- `nls_confint_normal`

These are numeric array routines, not S3 methods.
