# API map

## Direct numerical counterparts

| Upstream PMwR interface | Fortran interface |
|---|---|
| `.returns`, `returns.default` without dates | `simple_returns` |
| periodic returns | `period_returns` using integer period IDs |
| `returns_rebalance` | `portfolio_returns_weights` |
| `returns_position` | `portfolio_returns_positions` |
| `rebalance` | `rebalance_portfolio` |
| `valuation.position`, `pv` | `value_positions` |
| `unit_prices` | `unit_prices` |
| `div_adjust` | `dividend_adjust` |
| `split_adjust` | `split_adjust` |
| `journal` | `make_journal` and `journal_type` |
| journal sorting and concatenation | `sort_journal`, `append_journal` |
| `position.journal` | `positions_at` |
| `valuation.journal`, `jcf` | `journal_cashflows` |
| `.pl` | `pl_summary` |
| `.pl_stats` | `pl_path` |
| `split_trades` | `split_trade_runs` |
| `scale_to_unity` | `scale_trades_to_unity` |
| `close_on_first` | `close_on_first` |
| `limit` | `limit_position` |
| `tw_exposure` | `time_weighted_exposure` |
| `drawdowns` | `compute_drawdowns` |
| `streaks` | `compute_streaks` |
| `scale1` core scaling | `scale_to_level` |
| `rc(..., method="contribution")` | `link_return_contributions` |
| `rc(..., method="attribution")` | `return_attribution` |
| `is_valid_ISIN` | `valid_isin` |
| `is_valid_SEDOL` | `valid_sedol` |
| `quote32`, `q32` | `quote32_from_string`, `quote32_components`, `format_quote32` |

## Adapted interfaces

### `btest`

`run_backtest` receives an `n_time x n_asset` matrix of target positions or
weights. It supports signal masks, per-asset rebalance masks, execution lag,
partial rebalancing, transaction costs, cashflows, initial positions/cash, and
open-versus-close execution. It returns positions, suggested positions, cash,
wealth, cumulative costs, and a trade journal.

This preserves the accounting loop but does not execute arbitrary R callbacks
or maintain an R environment of strategy globals.

### `NAVseries`

The class is represented numerically by a level vector plus
`summarize_nav`. Date indexing and class metadata are left to the caller.

### `position`, `journal`, and account trees

The Fortran journal stores numeric timestamps and integer instrument/account
IDs. Hierarchical string-account expansion and R replacement methods are not
ported.

### `replace_weight`

`replace_group_weight` replaces a contiguous group and rescales a replacement
vector to the group's original total weight. R name-prefix behavior is omitted.

## Omitted non-computational or R-specific interfaces

- plotting and line methods
- HTML, Org, LaTeX, text, and sparkline reporting
- S3 print, summary formatting, subsetting, and coercion machinery
- `zoo`, `Date`, `POSIXct`, month/year formatting, and timestamp rounding
- parallel variations and replications
- browser links and packaged datasets
- arbitrary functions embedded in valuation and backtest calls
