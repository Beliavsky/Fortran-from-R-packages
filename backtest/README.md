# backtest-modern-fortran

A modern Fortran translation of the non-plotting computational surface of the
CRAN `backtest` package version 0.3-4.

The library implements portfolio-sort backtests using plain arrays. R data
frames, formulas, S4 classes, lattice plots, and `Date`/`factor` metadata are
replaced by explicit numeric matrices, integer group labels, and Fortran
derived types.

## Features

- R type-7 quantiles and numeric quantile bucketing
- Global or period-by-period signal buckets
- Numeric secondary-variable bucketing and pre-grouped categorical inputs
- Logical universe masks
- Multiple signal and return columns
- Grouped means, counts, trimmed means, and missing-value counts
- Return minimum, maximum, mean, median, sample standard deviation, and NA count
- Long-minus-short spreads and the original two-standard-error confidence band
- Natural-portfolio turnover, preserving the original no-price-change measure
- Overlapping holding-period weights and tri-bucket normalization
- Total and marginal counts
- Raw Sharpe ratios
- Cumulative quantile returns and worst-drawdown calculations

## Minimal example

```fortran
use backtest_kinds, only : dp
use backtest_engine, only : backtest_config, backtest_result, run_numeric_backtest

type(backtest_config) :: config
type(backtest_result) :: result
integer :: status

config%n_buckets = 5
config%by_period = .true.
config%natural = .true.

call run_numeric_backtest(signals, returns, config, result, status, &
                          period=period_codes, ids=security_ids)
```

`signals` and `returns` have observations in rows. `period_codes` and
`security_ids` are integer arrays aligned with those rows.

For already categorized signals, use `run_grouped_backtest` and pass bucket
codes directly.

## Command-line application

The CSV application expects:

```text
period,id,signal,return
```

Run:

```text
backtest_csv data/example_panel.csv 2 1
backtest_csv data/example_panel.csv 2 2
```

The second command uses two-period overlapping holdings.

## Building

```text
make check
make release-check
```

or with `fpm`:

```text
fpm test
fpm run demo_backtest
```

`fpm.toml` is provided, but validation for this release used the included
Makefile because `fpm` was not installed in the validation environment.

## Numerical conventions preserved

- Numeric buckets use R's default type-7 sample quantiles.
- Duplicate quantile breaks are reported as an error instead of silently
  assigning rank buckets.
- Returns are multiplied by observation weights before trimming and grouping.
- Trimming uses strict bounds at the 0.25% and 99.75% quantiles by default.
- Counts include missing return observations; NA counts are reported separately.
- Turnover ignores price changes and divides total long/short weight changes by
  four, matching the original implementation.
- Overlap normalization retains the original positive, negative, and zero
  weight-group behavior.
- The confidence interval is the original package's heuristic:
  spread plus or minus two times the full-sample return standard deviation
  divided by the square root of the combined low/high count.

## Explicit exclusions

- S4 `backtest` objects and generic method dispatch
- R formulas and expression evaluation; use a logical universe mask instead
- Data-frame merging, factor labels, and date/time classes
- Plotting, lattice panels, fan plots, and graphical annotations
- The packaged `starmine` R dataset

All underlying arrays used by the plotting methods, including period spreads,
cumulative bucket returns, turnover, and drawdown statistics, are available
numerically.

## License

The original package declares `GPL (>= 2)`. This translation is licensed under
GPL-2.0-or-later. See `LICENSE` and the SPDX header in every Fortran file.
