# PortfolioTesteR for Modern Fortran

This project is a modern Fortran/FPM translation of the computational core of
[`PortfolioTesteR` 0.1.4](upstream/PortfolioTesteR-main), an educational toolkit
for portfolio signals, construction, backtesting, performance analysis, machine
learning, parameter search, and walk-forward validation.

The port uses explicit numeric arrays and derived types rather than R data
frames, `xts`/`zoo` objects, S3 methods, or dynamically evaluated functions.
Unless a routine states otherwise, matrices are arranged as **time x asset**.
Feature panels are **time x asset x feature**, and parameter-grid combinations
are stored by column.

## Implemented computational areas

- Price momentum, distance, moving averages, RSI, stochastic oscillators, CCI,
  StochRSI, Bollinger bands, ATR, rolling correlation, and four rolling
  volatility estimators.
- Cross-sectional selection by rank, threshold, range, or percentile, with
  logical combinations and regime masks.
- Equal, signal, rank, inverse-volatility, hierarchical risk parity, equal-risk
  contribution, and maximum-diversification weights.
- Exposure, group, and turnover controls.
- Drift-aware share and cash accounting, integer or fractional shares,
  transaction costs, stop-loss exits, portfolio values, returns, and turnover.
- Return, volatility, Sharpe, Sortino, Omega, VaR, CVaR, drawdown, Calmar,
  benchmark, recovery-time, and regime statistics.
- Cross-sectional ranks, spreads, group-relative signals, breadth, and rolling
  correlation dispersion.
- Panel lags, interactions, labels, score transformations, pooled linear/ridge
  models, rolling predictions, information coefficients, bucket returns, and
  coverage diagnostics.
- Parameter-grid evaluation and walk-forward optimization through typed Fortran
  procedure callbacks.

See [`API_MAP.md`](API_MAP.md) for detailed coverage and
[`PORTING_NOTES.md`](PORTING_NOTES.md) for semantic differences.

## Build with FPM

```text
fpm build
fpm test
fpm run
```

The default executable is `demo_portfolio_tester`.

## Build with GNU Fortran scripts

Unix-like systems:

```text
scripts/build_and_test.sh checked
scripts/build_and_test.sh release
```

Windows command prompt:

```text
scripts\build_and_test.bat checked
scripts\build_and_test.bat release
```

The scripts compile and run all regression programs, examples, and the demo.

## Minimal example

```fortran
program momentum_example
  use portfolio_tester
  implicit none
  real(dp), allocatable :: prices(:,:), signal(:,:), selected(:,:), weights(:,:)
  type(backtest_result) :: result

  call generate_sample_prices(260, 12, prices, 123_i8)
  call calc_momentum(prices, 20, signal)
  call filter_top_n(signal, 4, selected)
  call weight_equally(selected, weights)
  call run_backtest(prices, weights, 100000.0_dp, result, &
    cost_bps=5.0_dp, integer_shares=.false., frequency=52.0_dp)

  print '(a,f10.4)', 'Total return: ', result%total_return
  print '(a,f10.4)', 'Sharpe:       ', result%sharpe
end program momentum_example
```

## License and provenance

The upstream package and this translation are distributed under the MIT
License. The complete upstream source snapshot is retained in `upstream/`.
See [`LICENSE`](LICENSE) and [`NOTICE.md`](NOTICE.md).

This software is for research and education. It is not investment advice.
