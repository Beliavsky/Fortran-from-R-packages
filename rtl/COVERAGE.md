# Computational coverage

## Translated

| Upstream routine | Fortran API | Notes |
|---|---|---|
| `GBSOption` and C++ `gbs` | `gbs_option`, `GBSOption` | Price and six Greeks; generalized cost of carry. |
| `CRROption` | `crr_option`, `CRROption` | European and American calls and puts. |
| `CRReuro` | `crr_euro_tree`, `CRReuro` | Returns full asset and option trees. |
| `spreadOption` | `spread_option`, `spreadOption` | Kirk-style price and Greeks. |
| `barrierSpreadOption` | `barrier_spread_option`, `barrierSpreadOption` | Up-and-out calls and down-and-out puts. |
| `simGBM` | `sim_gbm`, `simGBM` | Optional supplied increments and deterministic seed. |
| `simOU` | `sim_ou`, `simOU` | Euler process with optional supplied innovations. |
| `simOUt` | `sim_ou_time`, `simOUt` | Linearly interpolated time-varying mean. |
| `simOUJ` | `sim_ou_jump`, `simOUJ` | Compound-Poisson lognormal jumps. |
| `fitOU` | `fit_ou`, `fitOU` | Upstream closed-form OLS estimator. |
| `simMultivariates` | `sim_multivariates`, `simMultivariates` | Absolute changes, Kendall correlation, PSD safeguard, Gaussian simulation. |
| `bond` | `bond_value`, `bond` | Price, cash-flow arrays, and Macaulay duration. |
| `npv` | `npv_value`, `npv` | Discount-curve interpolation and flat-rate mode. |
| `swapIRS` | `interest_rate_swap`, `swapIRS` | Monthly, quarterly, semiannual, and annual schedules. |
| `swapFutWeight` | `swap_fut_weight`, `swapFutWeight` | Caller-supplied expiry date and holiday array. |
| `swapCOM` | `commodity_swap_prices`, `commodity_swap_from_calendar`, `swapCOM` | Direct weight or calendar-derived weight. |
| `efficientFrontier` | `efficient_frontier`, `efficientFrontier` | Random nonnegative fully invested portfolios. |
| `refineryLP` | `refinery_lp`, `refineryLP` | Native simplex solver replaces `lpSolve`. |
| `promptBeta` | `prompt_beta`, `promptBeta` | Full, bull, and bear betas. |
| `tradeStats` | `trade_stats`, `tradeStats` | Return, risk, Sharpe, Omega, participation, and drawdowns. |
| `returns` | `compute_returns`, `returns` | Absolute, relative, and log returns on arrays. |
| `rolladjust` | `roll_adjust_mask`, `rolladjust` | Identifies observations immediately after expiry for removal. |
| `tradeStrategySMA` | `moving_average_strategy`, `trade_strategy_sma`, `tradeStrategySMA` | Signals, trades, positions, returns, and cumulative equity. |
| `tradeStrategyDY` | `trade_strategy_dy`, `tradeStrategyDY` | Compatibility alias; upstream is identical to SMA strategy. |
| Numerical part of `chart_PerfSummary` | `trade_stats` | Plotting omitted. |

## Explicitly excluded

- `garch` and all `rugarch` integration.
- ZEMA, EIA, Genscape, Bank of Canada, and GIS downloads.
- Authentication, HTTP, JSON, and web API code.
- `swapInfo`, which is primarily API retrieval and plotting.
- `eia2tidy`, `eia2tidy_all`, and other ingestion wrappers.
- Plotting routines and Plotly/ggplot objects.
- `chart_zscore` STL decomposition, delegated upstream to `feasts`.
- R `xts`, `zoo`, `tsibble`, tibble, tidyverse, Rcpp registration, and serialization infrastructure.
- Bundled `.rda`, spreadsheet, Feather, RDS, and market-data files as compiled inputs.

The unmodified supplied package remains under `original/RTL-1.3.9` for provenance and comparison.
