# API

## Kinds and statuses

- `dp`: double-precision real kind.
- `int64`: 64-bit integer kind used for seeds.
- `jdmbs_success`
- `jdmbs_invalid_argument`
- `jdmbs_nonfinite_input`
- `jdmbs_numerical_warning`

## `type(jdmbs_control)`

- `day`: number of simulated calendar days.
- `monte_carlo`: number of paths per asset.
- `seed`: deterministic portable RNG seed.
- `days_per_year`: defaults to 365.
- `discount_rate`: continuous annual discount rate; defaults to zero.
- `legacy_mode`: defaults to true for upstream formula compatibility.
- `store_paths`: allocate and return `paths` when true.

## `type(jdmbs_result)`

- `call_price(:)`, `put_price(:)`: estimated option values.
- `call_se(:)`, `put_se(:)`: Monte Carlo standard errors.
- `terminal_price(:, :)`: asset by path terminal prices.
- `paths(:, :, 0:day)`: optional stored paths.
- `jump_events`: generated jump events. Shared events are counted once.
- `clipped_exponentials`: exponential evaluations clipped to avoid overflow.
- `status`, `message`: result status.

## `normal_bs`

```fortran
call normal_bs(start_price, mu, sigma, strike, result, control)
```

All input vectors must have the same size.

## `jdm_bs`

```fortran
call jdm_bs(start_price, mu, sigma, lambda, strike, result, control)
```

`lambda` is the expected number of jumps over the simulated contract horizon,
matching the effective interpretation of the upstream rescaling.
Each asset/path receives an independent Poisson jump stream. Jump multipliers
are uniform on `[0, 2]`, as in the R code.

## `jdm_new_bs`

```fortran
call jdm_new_bs(correlation_matrix, start_price, mu, sigma, lambda, &
   strike, result, control)
```

The matrix is retained under the upstream argument name, but numerically it is
a **jump-transmission matrix**, not a Brownian correlation matrix. For an event
originating at company `s`, asset `i` is multiplied by

```text
1 + jump_size * correlation_matrix(s, i)
```

where `jump_size` is uniform on `[-1, 1]`. Entries must lie in `[-1, 1]`.
