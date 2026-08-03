# PMwR-fortran

A modern Fortran/FPM translation of the computational core of the R package
`PMwR` 1.2-0 (Portfolio Management with R).

The port focuses on portfolio arithmetic and analysis rather than R's S3 class,
reporting, plotting, and date-time infrastructure. Public interfaces use
explicit vectors, matrices, logical masks, integer instrument identifiers, and
small derived types.

## Implemented functionality

- Simple and lagged returns for vectors and matrices
- Returns by integer period identifiers
- Portfolio returns from weights or positions
- Holdings and per-asset return contributions
- Portfolio rebalancing with weights, notionals, multipliers, truncation, and partial execution
- Position valuation and unit-price accounting for external cashflows
- Dividend and split adjustment
- Journal creation, sorting, concatenation, cashflow conversion, and position reconstruction
- Closed-trade P/L and average-cost realized/unrealized P/L paths
- Trade splitting at position sign changes, trade scaling, first-close logic, position limits, and time-weighted exposure
- Drawdown episodes, market streaks, and NAV summaries
- Geometric and logarithmic linking of return contributions
- Brinson, top-down, and bottom-up return attribution
- Matrix-driven portfolio backtesting with signal lag, per-asset rebalance masks, transaction costs, cashflows, open/close execution, and weight conversion
- ISIN and SEDOL check-digit validation
- US Treasury 32nds quote parsing and formatting

## Build with FPM

```text
fpm build
fpm test
fpm run demo_pmwr
fpm run --example backtest_targets
```

The manifest uses the FPM-compatible numeric version `1.2.0`, corresponding to
upstream version `1.2-0`.

## Build without FPM

On Unix-like systems:

```text
./scripts/build_all.sh check
./scripts/build_all.sh release
```

On Windows with GNU Fortran:

```text
scripts\build_all.bat check
scripts\build_all.bat release
```

## Minimal example

```fortran
program example
   use pmwr, only : dp, simple_returns
   implicit none
   real(dp), allocatable :: r(:)

   call simple_returns([100.0_dp, 105.0_dp, 102.0_dp], r)
   print *, r
end program example
```

## Backtesting interface

The upstream `btest()` evaluates user-supplied R callbacks at every timestamp.
`run_backtest` instead accepts a complete target-position or target-weight
matrix. This is deterministic, compiler-friendly, and suitable for strategy
signals computed by another Fortran routine. A lag shifts targets before
execution, and logical masks control signal updates and per-asset rebalancing.

## Scope notes

The following R-specific facilities are not reproduced:

- S3 classes and method dispatch for `NAVseries`, `journal`, `position`,
  `pricetable`, `pl`, and `btest`
- `zoo`, `Date`, and `POSIXct` processing
- formula/data-frame conveniences and instrument/account names
- plotting, HTML/Org/LaTeX/text reporting, sparklines, and browser integration
- parallel variation and replication orchestration
- arbitrary R callback execution inside backtests
- packaged R datasets

See `docs/API_MAP.md` and `docs/PORTING_NOTES.md` for detailed mappings and
adaptations.

## License

GPL-3.0-only. See `LICENSE`, `COPYING`, and `NOTICE.md`.
