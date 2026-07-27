# BCC1997 modern Fortran port

A self-contained modern Fortran and FPM translation of the computational
content of the R package `BCC1997` 0.1.1.

The library prices European calls and puts using the universal transform
solution of Bakshi, Cao, and Chen (1997). The model combines stochastic
variance, stochastic short rates, correlation between the asset and variance,
and compound-Poisson lognormal jumps.

## Build

```text
fpm build
fpm test
fpm run bcc1997_demo
fpm run --example basic_pricing
fpm run --example strike_curve
```

The project has no external numerical dependencies.

## Main interface

The familiar scalar interface is retained:

```fortran
use bcc1997
implicit none

type(bcc_result) :: result

result = bcc(kappav=1.5_dp, kappar=0.4_dp, thetav=0.04_dp, &
   thetar=0.03_dp, sigmav=0.30_dp, sigmar=0.10_dp, &
   muj=-0.05_dp, sigmaj=0.20_dp, rho=-0.60_dp, lambda=0.20_dp, &
   s0=100.0_dp, k=105.0_dp, v0=0.04_dp, r0=0.025_dp, t=1.25_dp)

print *, result%call, result%put
```

For reusable parameter sets, use the typed interface:

```fortran
type(bcc_parameters) :: parameters
type(bcc_result) :: result

parameters = bcc_parameters(spot=100.0_dp, strike=100.0_dp, &
   variance0=0.04_dp, rate0=0.01_dp, maturity=1.0_dp)
result = bcc_price(parameters)
```

`bcc_price_strikes` prices a vector of strikes. `integration_settings` controls
absolute and relative tolerances, panel width, maximum transform bound, and
adaptive depth.

## Result diagnostics

In addition to call and put values, `bcc_result` returns:

- the two transform probabilities;
- raw Fourier integrals and estimated errors;
- final integration bounds and function-evaluation counts;
- convergence status and a diagnostic message.

## Numerical method

The original R code calls `stats::integrate` on two improper Fourier
integrals. This port uses adaptive 15-point Gauss-Kronrod integration on
successive finite panels and terminates when several consecutive tail panels
are below the requested tolerance. Complex square roots and logarithms use
the Fortran principal branches, matching R's default complex arithmetic.

Stable series evaluation is used for `1 - exp(-z)` near zero. Terms multiplied
by an exactly zero long-run level or jump intensity are skipped before
potentially singular intermediate expressions are formed. This is important
for the two examples supplied with the original package.

## Parameter restrictions

The implementation follows the original package's note that variance and rate
volatilities must be positive. It also validates positive spot and strike,
nonnegative maturity, variance, mean-reversion speeds, long-run levels and jump
intensity, `rho` in `[-1, 1]`, and `1 + mu_j > 0`.

## Scope

The upstream package exports one computational routine, `BCC`; it is fully
translated. R documentation and list presentation are represented by typed
Fortran APIs. There is no plotting or data-download functionality in the
original package.

## License and citation

The upstream package declares `GPL (>= 2)`. This port is therefore licensed as
`GPL-2.0-or-later`. See `LICENSE`, `NOTICE`, and the unmodified source under
`original/BCC1997-0.1.1`.

Users of the model should cite Bakshi, Cao, and Chen (1997), DOI
`10.1111/j.1540-6261.1997.tb02749.x`, as requested by the original package.
