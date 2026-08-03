# API

All public procedures are available through:

```fortran
use greeks
```

All real arguments use `dp = kind(1.0d0)`.

## Result type

`type(greek_result)` contains:

- `status`: zero on success.
- `message`: diagnostic text on failure.
- `iterations`: iteration count where relevant.
- `names(:)`: requested Greek names.
- `values(:)`: estimates in the same order.
- `standard_errors(:)`: Monte Carlo standard errors, zero for exact methods.

## Exact European options

```fortran
call bs_european_greeks(spot, strike, rate, time, sigma, dividend, &
  payoff, requested, result)
```

Payoffs:

- `call`, `put`
- `cash_or_nothing_call`, `cash_or_nothing_put`
- `asset_or_nothing_call`, `asset_or_nothing_put`
- `digital_call`, `digital_put` as aliases for cash-or-nothing

Greeks:

`fair_value`, `delta`, `vega`, `theta`, `rho`, `epsilon`, `lambda`,
`gamma`, `vanna`, `charm`, `vomma`, `veta`, `vera`, `speed`, `zomma`,
`color`, and `ultima`.

`bs_european_price(...)` returns only the price.

## Exact geometric Asian options

```fortran
call bs_geometric_asian_greeks(spot, strike, rate, time, sigma, dividend, &
  payoff, requested, result)
```

Supports `call` and `put`, and the requested names `fair_value`, `delta`,
`rho`, `vega`, `theta`, `gamma`, and `vomma`.

## American options

```fortran
call binomial_american_greeks(spot, strike, rate, time, sigma, dividend, &
  payoff, requested, result, steps, eps)
```

The tree value is corrected by adding the exact European value and subtracting
the tree's European value, matching the original R wrapper. Supported Greeks
are `fair_value`, `delta`, `gamma`, `vega`, `theta`, `rho`, and `epsilon`.

`binomial_values(...)` exposes the raw American and European tree values.

## Malliavin European options

```fortran
call malliavin_european_greeks(spot, strike, rate, time, sigma, payoff, &
  requested, result, paths, seed, antithetic, payoff_fn)
```

Supports `fair_value`, `delta`, `vega`, `theta`, `rho`, and `gamma`.
The original routine has no dividend-yield argument. `payoff_fn` is an optional
pure callback with signature `(x, strike) result(value)`.

## Malliavin arithmetic Asian options

```fortran
call malliavin_asian_greeks(spot, strike, rate, time, sigma, dividend, &
  payoff, requested, result, model, jump_intensity, jump_scale, steps, &
  paths, seed, antithetic, payoff_fn, derivative_fn)
```

Models are `black_scholes` and `jump_diffusion`. The built-in jump-size law is
Student t with 3 degrees of freedom, multiplied by `jump_scale`, matching the
upstream default callback.

Supported names are `fair_value`, `delta`, `delta_d`, `rho`, `rho_d`,
`theta`, `theta_d`, `vega`, `vega_d`, `gamma`, and `gamma_kombi`.
Derivative-based estimators with a custom payoff require `derivative_fn`.

## Malliavin geometric Asian options

```fortran
call malliavin_geometric_asian_greeks(spot, strike, rate, time, sigma, &
  dividend, payoff, requested, result, model, jump_intensity, jump_scale, &
  steps, paths, seed, antithetic, payoff_fn)
```

Supports `fair_value`, `delta`, `rho`, `vega`, `theta`, and `gamma`.

## Control-variate arithmetic Asian options

```fortran
call bs_malliavin_asian_greeks(spot, strike, rate, time, sigma, dividend, &
  payoff, requested, result, steps, paths, seed)
```

Supports `call` and `put`, and `fair_value`, `delta`, `rho`, and `vega`.
It uses the exact geometric-Asian estimator as a control variate.

## Implied volatility

```fortran
call bs_implied_volatility(option_price, spot, strike, rate, time, dividend, &
  payoff, volatility, status, iterations, start_volatility, precision, max_iter)
```

uses Halley's method for European calls and puts.

```fortran
call implied_volatility(option_price, spot, strike, rate, time, dividend, &
  model, option_type, payoff, volatility, status, iterations, &
  start_volatility, precision, max_iter, steps, paths, seed)
```

uses the exact European solver when possible and a Newton update through the
general dispatcher otherwise.

## Dispatcher

```fortran
call option_greeks(spot, strike, rate, time, sigma, dividend, model, &
  option_type, payoff, requested, result, steps, paths, seed, antithetic)
```

Option types are `European`, `Digital`, `American`, `Asian`, and
`Geometric Asian`; matching is ASCII case-insensitive.

## Numerical helpers

The translated internal helpers are public for testing and advanced use:

- `make_bm`, `row_cumsums`
- `calc_i`, `calc_i_1`, `calc_i_2`, `calc_i_3`
- `calc_x`, `calc_log_x`, `calc_xw`, `calc_txw`
