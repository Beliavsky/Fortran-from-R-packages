# Dowd for modern Fortran

A dependency-free modern Fortran/FPM translation of the computational core of
R package **Dowd 0.12**, which contains market-risk examples based on Kevin
Dowd's *Measuring Market Risk*.

The project uses Fortran 2018, `implicit none`, a single double-precision kind,
standard allocatable arrays, and no external numerical libraries.

## Included numerical areas

- Normal, Student-t, lognormal, and log-Student-t VaR and expected shortfall
- Historical simulation, kernel-density, bootstrap, Box-Cox, and
  Cornish-Fisher risk measures
- Gumbel, Frechet, generalized Pareto, Hill, and Pickands tail methods
- Variance-covariance, adjusted variance-covariance, incremental hotspot, and
  PCA portfolio risk
- Binomial, Christoffersen, Lopez, Blanco-Ihle, Jarque-Bera, Kolmogorov-Smirnov,
  Kuiper, and Anderson-Darling backtesting/statistical routines
- Black-Scholes option prices and option VaR/ES simulation
- American put pricing and probability-weighted binomial VaR/ES
- Product, Gaussian, and Gumbel copulas and two-asset sum/VaR calculations
- Insurance, stop-loss, filter-strategy, defaultable-bond, and pension examples
- Normal and Student-t probability functions, quantiles, special functions,
  random generation, sorting, moments, and quantiles needed by the library

Plot-only R wrappers, R graphics, R formula/class dispatch, and package-specific
console presentation are intentionally omitted. See `PORTING.md` for the exact
scope and corrections made to inherited formulas.

## Build with FPM

```text
fpm build
fpm test
fpm run
fpm run --example option_risk
```

The default executable is `dowd_demo`.

## Build directly with GNU Fortran

```text
mkdir build
gfortran -std=f2018 -J build -I build -c \
  src/dowd_kinds.f90 src/dowd_math.f90 src/dowd_risk.f90 \
  src/dowd_portfolio.f90 src/dowd_backtests.f90 src/dowd_options.f90 \
  src/dowd_copulas.f90 src/dowd_simulations.f90 src/dowd.f90
ar rcs build/libdowd.a *.o
gfortran -std=f2018 -J build -I build test/test_dowd.f90 \
  build/libdowd.a -o build/test_dowd
./build/test_dowd
```

## Minimal use

```fortran
program example
  use dowd, only: dp, normal_var, normal_es
  implicit none

  print *, normal_var(0.001_dp, 0.02_dp, 0.99_dp, 10.0_dp)
  print *, normal_es (0.001_dp, 0.02_dp, 0.99_dp, 10.0_dp)
end program example
```

The library treats input series to historical routines as **profit/loss**:
positive values are profits and negative values are losses. Reported VaR and ES
are positive loss amounts when the lower tail contains losses.

## Modules

- `dowd_kinds`: precision constants
- `dowd_math`: distributions, special functions, statistics, RNG, sorting
- `dowd_risk`: univariate, nonparametric, bootstrap, and EVT risk measures
- `dowd_portfolio`: portfolio, hotspots, and PCA calculations
- `dowd_backtests`: coverage, independence, score, and goodness-of-fit tests
- `dowd_options`: European and American option risk
- `dowd_copulas`: bivariate copulas and sum distributions
- `dowd_simulations`: insurance, trading-rule, bond, and pension simulations
- `dowd`: umbrella module exporting the public library

## License

The original package says `License: GPL` without a version-specific notice. The
translation is offered under **GPL-2.0-only OR GPL-3.0-only**. See `LICENSE`,
`LICENSES/`, and `NOTICE.md`.
