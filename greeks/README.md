# greeks-fortran

A dependency-free modern Fortran/FPM translation of the computational layer of
R package `greeks` 1.5.6.

The library calculates analytical and Monte Carlo sensitivities for European,
American, arithmetic Asian, and geometric Asian options. It also provides
safeguarded implied-volatility solvers and callback interfaces for custom
payoffs and jump distributions.

## Build

```text
fpm build
fpm test
fpm run greeks_demo
fpm run --example analytic_greeks
fpm run --example malliavin_greeks
```

GNU Fortran users can run the supplied direct validation scripts:

```text
scripts/validate.sh
scripts\validate.bat
```

## Main API

```fortran
use greeks

type(greeks_result) :: result

result = bs_european_greeks(
  spot, strike, rate, time, volatility, dividend_yield, payoff_call)
```

`greeks_result` includes:

- `fair_value`, `delta`, `vega`, `theta`, `rho`, and `epsilon`
- `elasticity`, corresponding to the upstream Greek named `lambda`
- `gamma`, `vanna`, `charm`, `vomma`, `veta`, and `vera`
- `speed`, `zomma`, `color`, and `ultima`
- direct-payoff Asian estimators `delta_d`, `rho_d`, `theta_d`, `vega_d`, and
  `gamma_kombi`
- `ok` and `message` status fields

Payoff constants are:

```fortran
payoff_call
payoff_put
payoff_cash_call
payoff_cash_put
payoff_asset_call
payoff_asset_put
```

### Analytical European options

```fortran
result = bs_european_greeks(
  100.0_dp, 105.0_dp, 0.03_dp, 1.5_dp, 0.25_dp, 0.01_dp,
  payoff_call)
```

This translates all upstream European and binary formulas, including the
higher-order Greeks.

### Geometric Asian options

```fortran
result = bs_geometric_asian_greeks(
  100.0_dp, 100.0_dp, 0.02_dp, 1.0_dp, 0.25_dp, 0.0_dp,
  payoff_call)
```

The continuous geometric-average formulas use the risk-free discount rate and
match the upstream implementation.

### American options

```fortran
result = binomial_american_greeks(
  100.0_dp, 100.0_dp, 0.05_dp, 1.0_dp, 0.20_dp, 0.0_dp,
  payoff_put, steps=500)
```

The price uses a Cox-Ross-Rubinstein tree with the upstream European-tree bias
correction. Greeks other than value are finite differences around this corrected
price.

### Malliavin Monte Carlo

```fortran
type(mc_greeks_result) :: mc

mc = malliavin_european_greeks(
  100.0_dp, 100.0_dp, 0.03_dp, 1.0_dp, 0.20_dp,
  payoff_call, paths=40000, seed=123, antithetic=.true.)
```

`mc%estimate` contains the estimates and `mc%standard_error` contains pathwise
Monte Carlo standard errors.

Available routines are:

- `malliavin_european_greeks`
- `malliavin_geometric_asian_greeks`
- `malliavin_asian_greeks`
- `bs_malliavin_asian_greeks`, which uses the exact geometric Asian result as a
  regression control variate

Arithmetic and geometric Asian simulations support the package's uncompensated
compound-Poisson jump model. A custom jump generator can be supplied through a
`jump_sampler_callback`. The default jump law is Student t with three degrees
of freedom, matching the R default.

### Implied volatility

```fortran
type(implied_vol_result) :: iv

iv = bs_implied_volatility(
  option_price, spot, strike, rate, time, dividend_yield, payoff_call)
```

Additional routines are:

- `geometric_asian_implied_volatility`
- `american_implied_volatility`
- `implied_volatility`, accepting an arbitrary `price_callback`

The solvers combine bracketing with Halley or safeguarded Newton steps. The
result reports convergence, iteration count, and final pricing error.

## Design differences from R

- R named vectors and tibbles are typed Fortran result structures.
- A single Fortran call accepts scalar model parameters. Vectorized studies are
  expressed naturally with a caller loop or array operation.
- R function objects are typed procedure callbacks.
- R `NA` values and exceptions are replaced by explicit status fields.
- The random-number stream is deterministic but is not bitwise identical to
  `dqrng`.
- Shiny, plotly, ggplot2, and other user-interface code is excluded.

See `COVERAGE.md`, `PORTING_NOTES.md`, and `VALIDATION.md` for the exact mapping,
intentional changes, and test record.

## License

MIT, preserving the original package license and attribution. The unmodified R
package is under `original/greeks-1.5.6`.
