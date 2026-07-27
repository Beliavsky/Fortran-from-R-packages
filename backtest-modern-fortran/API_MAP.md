# API map

## Direct numerical translations

| Original R routine or method | Fortran procedure or field |
|---|---|
| `backtest()` numeric signals | `run_numeric_backtest` |
| `backtest.compute()` | `run_grouped_backtest` and `backtest_result` |
| `categorize()` numeric path | `categorize_quantiles` |
| period-by-period `categorize()` | `categorize_by_period` |
| categorical/factor signal path | precomputed bucket codes passed to `run_grouped_backtest` |
| `bucketize(..., mean)` | `bucketize_statistics` means |
| `bucketize(..., length)` | `bucketize_statistics` counts |
| NA bucketing | `bucketize_statistics` NA counts |
| `.bt.spread()` | `spread_statistics` |
| `.bt.mean()` | `mean_rows` |
| `.bt.sharpe()` | `sharpe_ratios` |
| `calc.turnover()` | `turnover_series` |
| `tribucket(..., scale=FALSE)` | `tribucket_weights` |
| `tribucket(..., scale=TRUE)` | `scale_period_weights` |
| `calc.true.weight()` | rolling accumulation inside `overlapping_weights` |
| `overlaps.compute()` | `overlapping_weights` plus `run_numeric_backtest` |
| `means()` | `backtest_result%means` |
| `counts()` | `backtest_result%counts` |
| `totalCounts()` | `total_counts` |
| `marginals()` | `marginal_counts` |
| `naCounts()` | `backtest_result%na_counts` |
| `turnover()` | `backtest_result%turnover` |
| `ci()` | `spread_statistics` confidence arrays |
| cumulative-return plot calculations | `cumulative_bucket_returns` |
| drawdown annotation calculations | `worst_drawdown` |
| return summary slot | `backtest_result%return_stats` |
| trimmed means slot | `backtest_result%trimmed_means` |

## Plain-array replacements

- R factor levels are integer bucket or group codes.
- R date values are integer period codes; original values are retained in
  `backtest_result%by_levels`.
- R universe expressions are logical masks passed as `universe=`.
- Multiple return and signal variable names are represented by column indices.
- S4 accessors are direct fields or procedures over `backtest_result`.

## Excluded infrastructure

- S4 class definitions, validity methods, generic dispatch, printing, and
  data-frame formatting
- Formulas and R expression evaluation
- Plotting and lattice/grid graphics
- `Date`, `POSIXt`, and factor-label metadata
- R `.RData` examples and the `starmine` dataset
