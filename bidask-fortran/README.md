# bidask-fortran

Modern Fortran/FPM translation of the numerical algorithms in `bidask` 2.1.5.
The library estimates bid-ask spreads from open, high, low, and close prices.
It is self-contained and has no external numerical dependencies.

## Implemented algorithms

- Efficient Discrete Generalized Estimator (`edge`)
- Fixed-width, adaptive-width, expanding, and endpoint EDGE estimates
- Abdi-Ranaldo (`AR`) and nonnegative two-period (`AR2`) estimators
- Corwin-Schultz (`CS`) and nonnegative two-period (`CS2`) estimators
- Roll estimator
- Generalized `OHL`, `OHLC`, `CHL`, and `CHLO` estimators
- Dot-separated averages such as `OHLC.CHLO`
- Multi-method `spread` dispatch
- Seeded OHLC simulation with observation probability, bid/ask bounce,
  intraperiod volatility, overnight volatility, and drift

R-specific `data.table`, `xts`, `zoo`, timestamps, and data-frame adapters are
not compiled. Numerical results are returned in typed arrays.

## Build

```text
fpm build
fpm test
fpm run bidask_demo
fpm run --example basic_estimators
fpm run --example rolling_and_simulation
```

## Basic use

```fortran
use bidask

real(dp) :: s
s = edge(open_prices, high_prices, low_prices, close_prices)
```

Multiple estimators can be evaluated through an `ohlc_data` object:

```fortran
type(ohlc_data) :: prices
type(spread_result) :: result
character(len=16) :: methods(3)

methods = ['EDGE            ', 'AR              ', 'OHLC.CHLO       ']
result = spread(prices, methods, na_rm=.true.)
```

The rolling result arrays have the same length as the input. Positions without
enough observations contain IEEE quiet NaNs.

## Simulation

```fortran
use iso_fortran_env, only: int64

type(ohlc_data) :: prices
prices = sim(1000, 390, spread=0.01_dp, volatility=0.03_dp, &
  seed=12345_int64)
```

`sim()` may also be called without `n` or `trades`; the upstream defaults of
10,000 periods and 390 trades per period are used.

## License and citation

MIT licensed. The upstream source, attribution, citation, and source archive
are retained in this project. See `LICENSE`, `NOTICE`, and `original/`.

Ardia, D., Guidotti, E., and Kroencke, T. A. (2024). Efficient estimation of
bid-ask spreads from open, high, low, and close prices. Journal of Financial
Economics, 161, 103916.
