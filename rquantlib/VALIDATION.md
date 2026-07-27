# Validation

## Environment

- GNU Fortran 14.2.0
- GNU ld/system linker
- LAPACK and BLAS from the Debian runtime environment
- Fortran standard mode: `-std=f2018`

## Debug flags

```text
-O0 -g -Wall -Wextra -Wimplicit-interface -Werror
-fcheck=all -fbacktrace
```

## Release flags

```text
-O2 -Wall -Wextra -Wimplicit-interface -Werror -fbacktrace
```

## Commands

```sh
./scripts/build_and_test.sh debug
./scripts/build_and_test.sh release
```

## Successful test output

```text
Option pricing and implied-volatility tests passed.
Date, curve, and bond tests passed.
SABR and Hull-White tests passed.
GPL-2.0-or-later source license checks passed.
debug build, tests, and applications passed.
release build, tests, and applications passed.
```

## Test coverage

The tests include:

- RQuantLib/QuantLib reference Black-Scholes prices and Greeks
- CRR American option reference pricing
- European cash-digital reference pricing
- Continuous geometric Asian reference pricing
- Barrier in/out parity
- Discrete-dividend behavior
- European and American implied-volatility recovery
- Arithmetic Asian Monte Carlo execution and error estimate
- Date round trips, weekdays, holidays, business-day conventions, schedules, and day counts
- Flat, explicit-zero, deposit/futures/swap bootstrap curve paths
- Zero and fixed-rate bond price/yield inversion
- Duration, convexity, cash-flow PV, floating coupons, caps, and floors
- Nelson-Siegel and Svensson synthetic curve recovery
- Off-ATM and ATM SABR formulas and synthetic calibration
- Hull-White time-zero discount-bond consistency
- Hull-White bond-option put-call parity
- Jamshidian swaption execution
- Hull-White caplet synthetic calibration
- All command-line applications and the standalone example

## Qualification

Passing these tests validates the implemented Fortran algorithms. It does not establish full QuantLib parity, historical named-calendar parity, or identical results from every engine and convention exposed by RQuantLib.
