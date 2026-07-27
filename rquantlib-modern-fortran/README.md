# rquantlib-modern-fortran

A standalone modern Fortran quantitative-finance library derived from the public numerical surface exposed by RQuantLib 0.4.28.

RQuantLib is primarily an R/Rcpp interface to the much larger C++ QuantLib library. This project does not embed or translate all of QuantLib. Instead, it provides independently implemented, plain-array Fortran procedures for the main self-contained pricing and term-structure workflows that could be validated without the external C++ library.

## Implemented numerical functionality

### Options

- Black-Scholes-Merton European calls and puts
- Value, delta, gamma, vega, theta, rate rho, and dividend rho
- Present-value adjustment for discrete cash dividends
- CRR binomial American calls and puts
- European and American cash-or-nothing and asset-or-nothing binaries
- Continuously monitored barrier approximations on a recombining tree
- Down/up, in/out barrier variants and in-out parity
- Continuous geometric-average Asian calls and puts
- Arithmetic-average Asian Monte Carlo with antithetic variates
- European and American implied-volatility inversion
- Strike-by-maturity European option arrays

### Dates, calendars, day counts, and schedules

- Gregorian date/serial-date conversion
- Day and month advancement, including end-of-month handling
- Weekend and user-supplied holiday calendars
- Following, modified-following, preceding, modified-preceding, and unadjusted conventions
- Business-day advancement and business-day counts
- Actual/360, Actual/365 Fixed, Actual/Actual ISDA, and 30/360 year fractions
- Forward or backward regular schedules

### Curves

- Flat continuously compounded curves
- Curves constructed from zero rates
- Log-linear interpolation and extrapolation of discount factors
- Discount factors, zero rates, and forward rates
- Simple deposit, futures, and par-swap bootstrapping
- Nelson-Siegel fitting by variable projection
- Nelson-Siegel-Svensson fitting by variable projection

### Bonds

- Zero-coupon price/yield conversion under simple, compounded, or continuous compounding
- Fixed-rate bond pricing from yield or a discount curve
- Fixed-rate yield inversion
- Cash-flow present value
- Macaulay and modified duration
- Convexity and one-basis-point sensitivity
- Floating-rate bond pricing from separate forward and discount curves
- Gearings, spreads, caps, floors, and redemption

### Volatility and short-rate models

- Hagan lognormal SABR implied volatility
- Fixed-beta SABR smile calibration
- Hull-White one-factor discount bonds
- Hull-White zero-coupon bond options
- Jamshidian European payer and receiver swaptions
- Hull-White caplets
- Two-parameter Hull-White calibration to caplet prices

## Build and test

GNU Fortran, LAPACK, and BLAS are required.

```sh
make debug
make release
```

The debug build uses runtime checking and backtraces. Both modes treat all compiler warnings as errors.

The project includes `fpm.toml`. A typical fpm build is:

```sh
fpm test --flag "-Wall -Wextra -Wimplicit-interface -Werror"
```

The supplied release was validated with the shell build harness because `fpm` was not installed in the validation environment.

## Applications

```sh
build/debug/bin/demo_rquantlib
build/debug/bin/price_option european call 100 100 0.01 0.04 1 0.2
build/debug/bin/price_option american put 100 105 0.01 0.04 1 0.25
build/debug/bin/price_bond 100 0.04 5 2 0.03
build/debug/bin/fit_curve data/example_yield_curve.csv ns
build/debug/bin/fit_curve data/example_yield_curve.csv svensson
```

## Important scope limits

The following RQuantLib functionality is not claimed:

- Exact QuantLib engine parity or use of QuantLib's C++ implementation
- Named national and exchange holiday calendars and their historical rule sets
- Full deposit/FRA/futures/swap/OIS helper and interpolation matrix
- Finite-difference, integral, Bjerksund-Stensland, Barone-Adesi-Whaley, and other alternate American engines
- Exact analytic barrier engines and all rebate timing conventions
- Discrete-fixing arithmetic Asian engines matching QuantLib
- Bermudan swaptions and the Black-Karasinski or G2++ engines
- Full SABR swaption cubes and market-convention calibration helpers
- Callable bonds
- Convertible zero, fixed-coupon, and floating-coupon bonds
- Full floating-rate coupon fixing, cap/floor volatility, and index-convention machinery
- Exponential-spline and polynomial fitted bond curves
- Hull-White calibration through complete cap/swap helper objects
- S3 classes, Rcpp modules, plotting, summaries, `zoo`, and R date metadata

Tree, Monte Carlo, calibration, and interpolation algorithms are tested numerical implementations, but exact QuantLib iteration paths, random streams, holiday decisions, and floating-point endpoints are not claimed.

## Licensing

The upstream RQuantLib package declares `GPL (>= 2)`. This project is distributed under **GPL-2.0-or-later**. Every Fortran source file contains a machine-checked SPDX identifier and GPL notice. `LICENSE` contains GNU GPL version 2.

The upstream QuantLib and Boost license texts are retained in `licenses/` for provenance. This standalone implementation does not link to either library.
