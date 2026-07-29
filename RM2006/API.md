# API

## Modules

### `rm2006_kinds`

Exports:

- `dp`: double-precision real kind, defined as `kind(1.0d0)`.

### `rm2006_module`

Exports the model procedures, defaults, and status codes.

## Main procedure

```fortran
call rm2006(data, covariance [, tau0] [, tau1] [, kmax] [, rho] [, status])
```

`rm2006` is a generic alias for `rm2006_covariance`.

### Arguments

- `data(:, :)`: input returns with shape `(n_obs, n_var)`.
- `covariance(:, :, :)`: allocatable output with shape
  `(n_var, n_var, n_obs+1)`.
- `tau0`: long time scale; default `1560.0_dp`; must exceed 1.
- `tau1`: shortest time scale; default `4.0_dp`; must be positive.
- `kmax`: number of time scales; default 14; must be at least 1.
- `rho`: geometric spacing ratio; default `sqrt(2.0_dp)`; must be positive.
- `status`: optional integer result code.

### Output indexing

- `covariance(:,:,1)` is the weighted multiscale backcast.
- `covariance(:,:,t+1)` is the estimate after return row `t`.
- `covariance(:,:,n_obs+1)` is the next-period forecast after all data.

The matrices are explicitly symmetrized at the end of the computation.

## Scale helper

```fortran
call rm2006_scale_weights(time_scales, weights [, tau0] [, tau1], &
                          [, kmax] [, rho] [, status])
```

Returns:

```text
time_scales(k) = tau1 * rho**(k-1)
raw_weight(k)  = 1 - log(time_scales(k)) / log(tau0)
weights        = raw_weight / sum(raw_weight)
```

This is the scale construction used by the original R function.

## Status codes

- `rm2006_success = 0`
- `rm2006_bad_shape = 1`
- `rm2006_bad_parameter = 2`
- `rm2006_nonfinite_data = 3`
- `rm2006_degenerate_weights = 4`

Use `rm2006_status_message(code)` to obtain a human-readable description.

On failure, the allocatable covariance output is returned with shape `(0,0,0)`.
