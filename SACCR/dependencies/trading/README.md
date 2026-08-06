# Trading modern Fortran port

This project translates the computational code of the R package `Trading` 3.2
into modern Fortran. It uses an FPM package layout, preserves the upstream
GPL-3 license and attribution, includes the upstream data files, and omits
plotting code.

The port covers:

- angular distance, Pearson correlation helpers, sample entropy, cross-sample
  entropy, normalized cross-sample entropy, variation of information,
  information-adjusted correlation, and information-adjusted beta;
- a native two-state Kalman filter, Rauch-Tung-Striebel smoother, and EM
  parameter estimator for dynamic beta and intercept;
- SA-CCR-style trade calculations, option supervisory deltas, CDO tranche
  deltas, maturity factors, supervisory durations, trade selection, FX dynamic
  setup, bond exposure splitting, CSA thresholds, collateral records, curves,
  and hash-table lookups;
- carbon footprint, carbon intensity, total carbon emissions, and weighted
  average carbon intensity;
- EuroMillions/EuroJackpot, Set For Life, UK Lotto, and Thunderball CSV loading,
  backtesting, payout/P&L calculations, top-number counts, and an iterator over
  all 139,838,160 Euro-lottery combinations;
- martingale run simulations and D'Alembert, Fibonacci, Labouchere, modified
  martingale, and specific-number roulette simulations;
- native CSV readers for the bundled curves, trades, CSA records, collateral,
  lottery histories, and example track record.

## Build with FPM

```text
fpm build
fpm test
fpm run trading_demo
```

The package has no external Fortran dependencies.

## Build with GNU Fortran

On Unix-like systems:

```text
./build_gfortran.sh
```

On Windows with `gfortran` available on `PATH`:

```text
build_gfortran.bat
```

## Minimal use

```fortran
program example
  use trading, only : dp, angular_distance, trade_t
  implicit none

  real(dp) :: returns(4, 2)
  real(dp) :: distance(2, 2)
  type(trade_t) :: swap

  returns(:, 1) = [0.01_dp, -0.02_dp, 0.015_dp, 0.005_dp]
  returns(:, 2) = [0.008_dp, -0.018_dp, 0.013_dp, 0.006_dp]
  call angular_distance(returns, distance)

  call swap%configure_class("IRDSwap")
  swap%notional = 10000.0_dp
  swap%si = 0.0_dp
  swap%ei = 10.0_dp
  swap%buy_sell = "Buy"

  write(*, '(f10.6)') distance(1, 2)
  write(*, '(f12.2)') swap%calc_adjusted_notional()
end program example
```

## Module layout

- `trading`: umbrella module
- `trading_dependence`: entropy and dependence metrics
- `trading_dynamic_beta`: Kalman filter, smoother, and EM
- `trading_trades`: trade representation and computational methods
- `trading_csa`: CSA and collateral types
- `trading_curve`: linear and natural-cubic interpolation
- `trading_hash_table`: two-column lookup tables
- `trading_climate`: portfolio climate metrics
- `trading_lottery`: loaders, backtests, P&L, and combination iterator
- `trading_betting`: stochastic betting strategy simulators
- `trading_io`: CSV readers
- `trading_stats`, `trading_strings`, `trading_kinds`: shared utilities

## Data

The upstream CSV files are copied into `data/`. Programs should run from the
package root when using paths such as `data/example_trades.csv`.

## Scope and compatibility

R data frames, reference classes, reflection, plotting, and package-attachment
messages do not have direct Fortran equivalents. They are represented by typed
arrays, derived types, type-bound procedures, CSV readers, and an umbrella
module. See `docs/API_MAP.md` and `docs/PORTING_NOTES.md` for exact mappings and
intentional corrections.

## License

The upstream package declares `GPL-3`. This translation is distributed under
`GPL-3.0-only`. The complete upstream R source and data used for the translation
are retained under `original-r/` for attribution and auditability.
