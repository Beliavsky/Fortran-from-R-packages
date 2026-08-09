# API

## Model callback

```fortran
subroutine model(x, par, y, ierr)
    real(dp), intent(in) :: x(:), par(:)
    real(dp), intent(out) :: y(:)
    integer, intent(out) :: ierr
end subroutine
```

Set `ierr=0` for a valid evaluation.  A nonzero value marks the evaluation as
invalid.

## `fit_onls`

```fortran
call fit_onls(model, x, y, start, result, control, lower, upper, weights, fixed)
```

Only the first five arguments are required.

- `x(:)`, `y(:)`: two-dimensional data.
- `start(:)`: parameter starting values.
- `lower`, `upper`: scalar or parameter-length bounds.
- `weights(:)`: nonnegative observation weights.
- `fixed(:)`: logical mask applied in the ONLS stage, matching upstream.

## `onls_control`

Important fields:

- `lm`: nested `lm_control`.
- `window = 12`.
- `extend = [0.2,0.2]`.
- `projection_tol = sqrt(epsilon(1.0_dp))`.
- `mimic_r_unsorted_weights = .true.`.

The nested `lm_control` contains `maxiter`, `maxfev`, `ftol`, `xtol`, `gtol`,
`fd_step`, and `lambda0`.

## `onls_result`

Key fields:

- `par_nls`, `par_onls`;
- sorted `x`, `y`, and the effective `weights`;
- nearest points `x0`, `y0`;
- `fitted_nls`, `fitted_onls`;
- `resid_nls`, `resid_onls`;
- `distance_o`: raw Euclidean point-to-curve distances;
- `resid_o`: upstream-style weighted orthogonal residuals;
- `rss_nls`, `rss_vertical`, `rss_orthogonal`;
- `gradient`, `covariance`, `stderr`;
- `ortho_angle`, `ortho`;
- solver iteration/status fields and `message`.

## Helpers

```fortran
llv = vertical_loglik(result)
llo = orthogonal_loglik(result)
r   = orthogonal_residuals(result)
```

`check_orthogonality` is public when a diagnostic must be recomputed explicitly.
