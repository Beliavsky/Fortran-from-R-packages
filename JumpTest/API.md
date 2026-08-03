# API reference

All public procedures are available from the aggregate module `jumptest`.
Real calculations use `dp = real64`.

## Result types

### `simulation_result`

- `price(:)`: simulated log-price path.
- `variance(:,:)`: latent volatility or factor paths.
- `jump_count(:)`: Poisson event counts when applicable.
- `jump_measure(:)`: upstream-compatible jump scale/measure.
- `jump_size(:)`: realized Gaussian jump increments.
- `seed`: normalized user seed supplied to the simulation.
- `status`: one of the `JT_*` status constants.

### `statp_result`

Scalar `stat`, `pvalue`, and `status` fields.

### `adjp_result`

Allocatable arrays `stat(:)`, `pvalue(:)`, and `adjp(:)`, plus `status`.

### `pcombine_result`

- `pvalue(:,:)`: periods by methods.
- `methods(:)`: normalized method labels.
- `status`.

## Simulation procedures

```fortran
call sv(intervals_per_period, periods, result, p0, mu, v0, b, alpha, sigma, seed)
call svj(intervals_per_period, periods, result, p0, lambda, mu, v0, b, &
         alpha, sigma, sigma1, seed)
call sv1f(intervals_per_period, periods, result, p0, mu, v0, beta0, beta1, &
          alphav, correlation, seed)
call sv1fj(intervals_per_period, periods, result, p0, lambda, mu, v0, beta0, &
           beta1, alphav, correlation, seed)
call sv2f(intervals_per_period, periods, result, p0, mu, v1, v2, beta0, &
          beta1, beta2, alpha1, alpha2, beta_v2, r1, r2, seed)
```

All model parameters after `result` are optional and use the upstream defaults.
Seeds are signed 64-bit integers.

Low-level deterministic recursions are also public:

```fortran
call lp_path(...)
call pvc_path(...)
call pv2_path(...)
```

These accept caller-supplied shocks and are useful for regression testing or
integration with another random-number generator.

## Jump statistics

```fortran
value = bns_statistic(returns, status)
value = amin_statistic(returns, status)
value = amed_statistic(returns, status)
call jumptestday(returns, result, method)
call jumptestperiod(return_matrix, result, method)
```

`return_matrix(:,j)` is one interval/day. Supported method labels are `BNS`,
`Amed`, and `Amin`, case-insensitively.

## Combining and pooling p-values

```fortran
call pcombine(return_matrix, methods, combined)
call ppool(pvalue_matrix, pooled, method)
```

`pcombine` requires at least two method names. `ppool` supports:

- `SD`: correlation-adjusted Stouffer.
- `FD`: correlation-adjusted Fisher.
- `SI`: independent Stouffer.
- `FI`: independent Fisher.
- `MI`: minimum p-value.
- `MA`: maximum p-value.

Input p-values are clipped to `[1e-5, 1-1e-5]`, matching the R package. Every
pooled result includes Benjamini-Hochberg adjusted p-values.

## Status constants

- `JT_SUCCESS`
- `JT_INVALID_ARGUMENT`
- `JT_INVALID_DIMENSION`
- `JT_NONFINITE_INPUT`
- `JT_DEGENERATE_SAMPLE`
- `JT_NUMERICAL_FAILURE`

`status_message(status)` returns a descriptive string.
