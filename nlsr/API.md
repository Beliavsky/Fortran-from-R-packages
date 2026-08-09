# API

## Main solver

### `nlfb(start, nres, resfn, result, ...)`

Arguments:

- `start(:)` - starting parameter vector.
- `nres` - residual-vector length.
- `resfn` - typed residual callback.
- `result` - `type(nlsr_result)` output.
- `control` - optional `type(nlsr_control)`.
- `jacfn` - optional analytic residual-Jacobian callback.
- `lower`, `upper` - scalar or parameter-length bounds.  Equal lower/upper values mask a parameter.
- `weights` - optional fixed nonnegative residual weights.
- `weightfn` - optional callback producing dynamic weights.

The solver minimizes the weighted sum of squared residuals.

## Control

`type(nlsr_control)` mirrors the upstream defaults where they are numerical:

- `femax = 10000`
- `jemax = 5000`
- `lamda = 1e-4`
- `laminc = 10`
- `lamdec = 4`
- `nbtlim = 6`
- `ndstep = 1e-7`
- `offset = 100`
- `phi = 1`
- `psi = 0`
- `rofftest = .true.`
- `smallsstest = .true.`
- `stepredn = 0`
- `scale_offset = 1`
- `jacobian_method = jac_central`

## Jacobian approximations

Constants:

- `jac_forward`
- `jac_backward`
- `jac_central`
- `jac_richardson`

Generic worker:

```fortran
call numerical_jacobian(resfn, par, res0, jac, method, ndstep, bdmask, feval, ierr)
```

Familiar wrappers matching upstream names are also exported:

- `jafwd`
- `jaback`
- `jacentral`
- `jand`

## Residual utilities

- `residual_gradient` - computes `2 J^T r` and residual sum of squares.
- `residual_sum_squares(r)` - returns `dot_product(r,r)`.

## Inference

```fortran
call nlsr_standard_errors(result, covariance, se, sigma, ok)
```

This uses the final weighted Jacobian and residual variance on unmasked parameters.

## Logistic self-start utilities

The computational parts of upstream `SSlogisJN` are exposed as:

- `logistic_value`
- `logistic_jacobian`
- `logistic_initial`

`logistic_initial` follows upstream's initialization: `Asym=2*max(y)`, followed by a linear regression of `log(Asym/y - 1)` on the input.

## Status constants

- `nlsr_ok`
- `nlsr_max_feval`
- `nlsr_max_jeval`
- `nlsr_no_progress`
- `nlsr_bad_input`
- `nlsr_callback_error`
- `nlsr_singular`
