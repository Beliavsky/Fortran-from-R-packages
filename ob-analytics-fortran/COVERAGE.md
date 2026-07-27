# Computational coverage

This file maps the R package `obAnalytics` 0.1.2 to the modern Fortran API.

## Fully translated numerical and data-processing routines

| R routine | Fortran routine | Notes |
|---|---|---|
| `vectorDiff` | `vector_diff` | First difference with a leading zero. |
| `reverseMatrix` | `reverse_matrix` | Reverses matrix rows. |
| `norml` | `normalize` | Optional caller-supplied limits. Constant inputs return zero. |
| `intervalSumBreaks` | `interval_sum_breaks` | Interval sums from cumulative break positions. |
| `vwap` | `weighted_average` | Volume-weighted mean. |
| `intervalVwap` | `interval_vwap` | Interval VWAP values. |
| `intervalPriceLevelGaps` | `interval_price_level_gaps` | Counts zero-volume levels by interval. |
| `loadEventData` | `read_event_csv` | Fixed seven-column CSV, sanitization, sorting, and fills. |
| `sMatrix` | `similarity_matrix_equal`, `similarity_matrix_time` | Typed similarity constructors. |
| `alignS` | `needleman_wunsch` | Global sequence alignment and matched-index output. |
| `eventMatch` | `event_match` | Nearest match followed by alignment when needed. |
| `matchTrades` | `match_trades` | Maker/taker, price, direction, and IDs. |
| `setOrderTypes` | `set_order_types` | All six original order categories. |
| `priceLevelVolume` | `price_level_volume` | Sparse price-level cumulative depth. |
| `filterDepth` | `filter_depth` | Opening-state clamp and closing zero records. |
| `depthMetrics` | `depth_metrics` | Best levels and configurable BPS bins. |
| `getSpread` | `get_spread` | Removes unchanged consecutive states. |
| `orderAggressiveness` | `order_aggressiveness` | BPS distance from preceding best quote. |
| `orderBook` | `reconstruct_order_book` | Active orders, cumulative liquidity, and BPS. |
| `tradeImpacts` | `trade_impacts` | Grouped market-order impact statistics. |
| `processData` | `process_data` | Complete pipeline, warm-up, and zombie filtering. |

## Typed replacements for R structures

R data frames and lists are represented by the derived types in `ob_types`.
Variable-length outputs use allocatable arrays. Factor levels are public integer
constants and have string conversion helpers.

## Excluded R-specific functionality

- `toZoo` and all `zoo`/time-series adapter behavior.
- `loadData` and `saveData`, which use RDS serialization.
- All `plot*` functions and `themeBlack`.
- `ggplot2`, `reshape2`, and R formula/list presentation infrastructure.
- Compiled use of `lob.data.RData`; the original data file is retained.

The numerical arrays required by the plotting routines are available from the
translated depth, spread, book, event, trade, and impact APIs.
