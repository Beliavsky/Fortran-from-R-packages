# API

All public procedures are available from `mcrp_module`.

## Kinds and status values

- `dp`
- `mcrp_success`
- `mcrp_invalid_shape`
- `mcrp_invalid_argument`
- `mcrp_numerical_failure`
- `mcrp_max_iterations`

## Centered co-moments

```fortran
m2(r [, status]) -> real(dp), allocatable :: matrix(:, :)
m3(r [, status]) -> real(dp), allocatable :: matrix(:, :)
m4(r [, status]) -> real(dp), allocatable :: matrix(:, :)
```

For `r(T,N)`:

- `m2` returns `N x N` and divides by `T-1`, matching R's covariance convention.
- `m3` returns `N x N^2` and divides by `T`.
- `m4` returns `N x N^3` and divides by `T`.

The flattened tensor ordering matches R's nested `kronecker` calls:

- `M3(i,(j-1)N+k) = E[x_i x_j x_k]`
- `M4(i,((j-1)N+k-1)N+l) = E[x_i x_j x_k x_l]`

## Raw portfolio moments

```fortran
pm2(r, w [, status]) -> scalar
pm3(r, w [, status]) -> scalar
pm4(r, w [, status]) -> scalar

dm2(r, w [, status]) -> vector(:)
dm3(r, w [, status]) -> vector(:)
dm4(r, w [, status]) -> vector(:)

cm2(r, w [, percentage, status]) -> vector(:)
cm3(r, w [, percentage, status]) -> vector(:)
cm4(r, w [, percentage, status]) -> vector(:)
```

`percentage` defaults to `.true.`. Percentage contributions sum to one;
absolute contributions sum to the corresponding raw moment.

## Standardized portfolio measures

```fortran
port_risk(r, w [, status])
port_risk_deriv(r, w [, status])
port_risk_contrib(r, w [, percentage, status])

port_skew(r, w [, status])
port_skew_deriv(r, w [, status])
port_skew_contrib(r, w [, percentage, status])

port_kurt(r, w [, status])
port_kurt_deriv(r, w [, status])
port_kurt_contrib(r, w [, percentage, status])
```

These correspond to the original `PortRisk`, `PortSkew`, and `PortKurt`
families. Fortran is case-insensitive, so the lowercase spellings naturally
serve as the original names.

## Multiple-criteria optimization

```fortran
subroutine mcrp(start, returns, result, lambda, active, lower, upper, &
   max_iterations, tolerance, initial_step)
```

Arguments:

- `start(N)`: raw starting parameters.
- `returns(T,N)`: return observations.
- `result`: `type(mcrp_result)`.
- `lambda(3)`: criterion weights for variance, skewness, and kurtosis.
- `active(3)`: explicitly enable or disable each criterion.
- `lower(N)`, `upper(N)`: optional bounds on raw parameters.
- `max_iterations`: default 4000.
- `tolerance`: default `1e-10`.
- `initial_step`: default `0.10`.

A NaN value in `lambda(i)` also disables criterion `i`, mirroring R's use of
`NA`. `active` is usually clearer in Fortran.

### `mcrp_result`

- `weights(:)`
- `raw_parameters(:)`
- `variance_contributions(:)`
- `skewness_contributions(:)`
- `kurtosis_contributions(:)`
- `objective`
- `iterations`
- `evaluations`
- `status`
- `converged`
- `message`

## Objective evaluation

```fortran
mcrp_objective_value(returns, x [, lambda, active, status]) -> scalar
```

This is useful for testing, diagnostics, and comparing candidate portfolios.
The objective is scale invariant in `x`.
