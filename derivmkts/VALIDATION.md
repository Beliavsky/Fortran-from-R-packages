# Validation report

## Environment

- GNU Fortran 14.2.0
- Fortran 2018 mode
- Linux x86-64
- No external numerical libraries

Strict debug flags:

```text
-std=f2018 -Wall -Wextra -Wconversion-extra -Wimplicit-interface
-Werror -fcheck=all -fbacktrace
```

An optimized `-O2` build was also compiled and run.

The translated project contains 2,156 lines across 27 Fortran source, test,
application, and example files.

## Test programs

```text
test_barriers_compound: PASS
test_binomial_jumps: PASS
test_black_scholes: PASS
test_greeks_misc: PASS
test_simulation_asian: PASS
```

The application and all examples also compile and run:

```text
derivmkts_demo
analytic_options
monte_carlo
simulation_and_tree
```

## Coverage of checks

- Fixed high-precision Black-Scholes call and put references
- Asset-binary and cash-binary decomposition identities
- Implied call and put volatilities
- Implied underlying prices
- Fixed bond price, yield, Macaulay duration, modified duration, and convexity
- Fixed geometric-average Asian price and strike references
- Every standard barrier in/out parity identity
- Every asset and cash binary barrier in/out parity identity
- Immediate and deferred rebate behavior
- Perpetual exercise-barrier sanity checks
- Bivariate normal value at the origin
- Fixed Geske compound-option prices and critical spots
- European binomial convergence to Black-Scholes
- American put dominance over its European counterpart
- Binomial terminal probability normalization
- Zero-intensity jump prices equal Black-Scholes prices
- Jump asset/cash decomposition and put-call parity
- Deterministic path reproduction from identical seeds
- Single- and multi-asset simulation dimensions and positivity
- Geometric Asian Monte Carlo agreement with analytical values
- Arithmetic Asian control-variate variance reduction
- Numerical Greeks against analytical Black-Scholes Greeks
- Quincunx count totals and theoretical binomial probabilities

## Portability safeguards

- No assumptions about short-circuit logical evaluation
- No fixed `kind=8`; `dp = kind(1.0d0)` is used throughout
- No external BLAS, LAPACK, R, `mnormt`, or plotting dependency
- Combined relative tolerances are used for platform-sensitive transcendental calculations
- All source files are free-form and remain within 132 columns
- All Fortran program units use `implicit none`

## FPM

FPM was not installed in the validation environment. The manifest was parsed as
TOML, and the project follows FPM's automatic `src`, `app`, `example`, and
`test` discovery conventions. Direct GNU Fortran scripts are included under
`scripts/` for independent validation on Linux and Windows.
