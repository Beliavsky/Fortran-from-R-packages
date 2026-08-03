# ragtop-fortran

A self-contained modern Fortran translation of the computational core of the
R package `ragtop` 2.0.0. The package prices equity derivatives and
credit-sensitive securities under Black-Scholes diffusion with optional jump
to default, including an equity-linked power-law hazard specification.

## Features

- Black-Scholes prices, delta, and vega with default and dividends
- Piecewise rate and volatility term structures
- Backward implicit finite-difference pricing
- European and American options
- Zero-coupon, coupon, callable, and convertible bonds
- Discrete dividends, coupons, calls, puts, recovery, and conversion
- Constant or stock-price-linked default intensity
- Implied-volatility and equivalent-volatility calibration
- Market intensity-link calibration
- Delta, gamma, vega, rho, hazard sensitivity, and theta

The implementation has no nonstandard library dependency.

## Build with FPM

```console
fpm build
fpm test
fpm run
```

Run an individual example with, for example:

```console
fpm run --example example_convertible_bond
```

## Basic use

```fortran
use ragtop

type(market_spec) :: market
real(dp) :: value

market%short_rate = 0.06_dp
market%volatility = 0.20_dp
value = american(put_option, 100.0_dp, 110.0_dp, 1.0_dp, &
                 market, min_steps=200)
```

For a convertible bond, create a `cashflow_schedule`, construct an
`instrument_spec` with `ConvertibleBond`, and call `find_present_value`.

## Array conventions

Fortran vectors use their natural one-dimensional representation. Price grids
store one value per stock node. Times are nonnegative year fractions. Rates,
volatilities, hazards, and dividend rates are in natural units, not percent.

## Licensing

The translation and retained upstream source are distributed under
GPL-2.0-or-later. See `LICENSE`, `NOTICE.md`, and `upstream/ragtop-master`.
