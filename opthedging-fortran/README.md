# OptHedging modern Fortran library

This project is a modern Fortran and FPM translation of `OptHedging` 1.0 by
Bruno Remillard. It computes option values and minimum-variance hedging
strategies under independently and identically distributed periodic returns.
The method follows Chapter 3 of *Statistical Methods for Financial Engineering*
by Bruno Remillard.

The original package is licensed under GPL version 2 or later. This derivative
work preserves that license as `GPL-2.0-or-later`.

## Implemented functionality

- Linear interpolation and extrapolation on a uniform price grid.
- Call and put payoffs on discounted prices.
- Backward Monte Carlo valuation under IID log excess returns.
- Period-by-period auxiliary hedge coefficients.
- Initial and subsequent optimal hedge positions.
- Typed results with status and diagnostic messages.
- Deterministic seeding and standard-normal simulation for examples.

The R `.C` interface and R list/matrix construction are not needed in Fortran.
All numerical functionality from the package is translated.

## Build with FPM

```text
fpm build
fpm test
fpm run opthedging_demo
fpm run --example basic_hedging
fpm run --example interpolation_and_rehedging
```

On Windows with GNU Fortran, the included direct validation script can also be
run from a Developer Command Prompt or ordinary `cmd.exe`:

```text
scripts\validate_gfortran.bat
```

## Basic use

```fortran
program example
   use opthedging, only : dp, hedging_iid, hedging_result
   implicit none

   real(dp) :: log_returns(6)
   type(hedging_result) :: fit

   log_returns = log([0.85_dp, 0.93_dp, 0.99_dp, 1.02_dp, 1.08_dp, 1.18_dp])
   fit = hedging_iid(log_returns, 0.5_dp, 100.0_dp, 0.03_dp, .true., &
      4, 401, 50.0_dp, 150.0_dp)
   if (.not. fit%ok) error stop fit%message

   print *, fit%option_value_at(1, 100.0_dp)
   print *, fit%initial_hedge_at(100.0_dp)
end program example
```

`log_returns` contains simulated periodic log excess returns. Internally the
simple excess returns are computed as `exp(log_returns) - 1`. Prices and option
values are discounted, matching the original package.

## Result fields

A `hedging_result` contains:

- `s`: uniform discounted-price grid.
- `c(period, grid_point)`: option values.
- `a(period, grid_point)`: auxiliary hedge coefficients.
- `phi1`: initial number of shares at each grid point.
- `rho`: return-moment coefficient used in the hedge.
- `discounted_strike`: `strike * exp(-rate * maturity)`.
- `ok` and `message`: validation status.

Type-bound methods evaluate the grid away from its nodes:

- `option_value_at(period, spot)`
- `auxiliary_at(period, spot)`
- `shares_at(period, spot, portfolio_value)`
- `initial_hedge_at(spot)`

The `period` index is one-based. Period 1 is time zero and period `n_periods`
is the last rebalancing date before maturity.

## Numerical method

For simple excess return `xi`, the method estimates

```text
rho = E[xi] / E[xi^2]
c2  = 1 - E[xi] * rho
z   = (1 - rho * xi) / c2
```

The terminal discounted payoff is propagated backward with the pricing kernel
`z`. At each date, the auxiliary coefficient is the sample mean of
`xi * next_value` divided by `E[xi^2]`. Values between or outside grid nodes use
linear interpolation or linear extrapolation, preserving the behavior of the
original package.

## Project layout

- `src`: Fortran library modules.
- `app`: demonstration program.
- `example`: focused examples.
- `test`: regression and mathematical-identity tests.
- `scripts`: direct GNU Fortran validation scripts.
- `original`: unmodified original package source.
- `provenance`: supplied archive and SHA-256 manifests.

See `COVERAGE.md`, `PORTING_NOTES.md`, and `VALIDATION.md` for details.
