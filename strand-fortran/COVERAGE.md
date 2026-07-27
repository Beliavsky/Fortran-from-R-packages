# Computational coverage

This document maps the numerical content of `strand` 0.2.3 to the Fortran port.
The R package is retained unmodified under `original/strand-0.2.3`.

## Utilities

| Upstream routine | Fortran implementation | Status |
|---|---|---|
| `.norm` | `rank_normal` | Implemented |
| `normalize` | `normalize`, `normalize_grouped` | Implemented for numeric arrays and integer group codes |
| `adjust` | `adjust_numeric` | Implemented for numeric regressors; categorical variables can be supplied as dummy columns |
| `drawdown` | `maximum_drawdown` | Implemented |

`rank_normal` follows the upstream sequence: standardize, round to 11 decimal
places for tie handling, average ranks, transform with the inverse normal CDF,
and rescale to the requested sample standard deviation.

## Cross-section and portfolio state

| Upstream class/method | Fortran implementation | Status |
|---|---|---|
| `CrossSection$update` | `update_cross_section` | Implemented |
| carry-forward values | `carry_columns`, `carry_values` arguments | Implemented |
| NA replacements | `replace_columns`, `replace_values` arguments | Implemented with IEEE NaN |
| `CrossSection$periodStats` | `compute_cross_section_stats` | Implemented |
| `Portfolio$initialize` | `initialize_portfolio` | Implemented |
| `Portfolio$getConsolidatedPositions` | `consolidated_shares` | Implemented |
| `Portfolio$applyAdjustmentRatio` | `apply_adjustment_ratio` | Implemented |
| internal/external positions | `portfolio_state` | Implemented |

Data-frame column maps and file-backed cross sections are R adapters and are not
compiled.

## Portfolio optimization

The upstream `PortOpt` R6 class is consolidated into `optimize_portfolio` and
typed configuration/result structures.

| Upstream numerical component | Fortran implementation | Status |
|---|---|---|
| alpha objective | `optimize_portfolio` | Implemented |
| full and half-way target-weight policies | `optimizer_config%target_weight_policy` | Implemented |
| explicit target weights | `strategy_spec%target_*_weight` | Implemented |
| maximum target-weight change | `optimizer_config%max_weight_change` | Implemented |
| strategy long/short market-value equalities | native LP constraints | Implemented |
| ADV/LMV/SMV position limits | `strategy_spec%position_limit_*` | Implemented |
| ADV trading limits | `strategy_spec%trading_limit_pct_adv` | Implemented |
| investability and no-side-switching rules | native variable bounds | Implemented |
| factor constraints | `factor_constraint` | Implemented |
| category constraints | `category_constraint` | Implemented |
| joint absolute-net-trade auxiliaries | native LP constraints | Implemented |
| turnover constraint | `optimizer_config%turnover_limit` | Implemented |
| constraint loosening | `optimizer_config%loosening_sequence` | Implemented |
| Rglpk/Rsymphony solve | `strand_simplex` | Replaced by native two-phase simplex |
| share rounding and post-round NMV | `optimization_result` | Implemented |
| max position/order diagnostics | `optimization_result%max_*` | Implemented |

The public R methods that expose sparse matrices, data frames, or R6 state are
not reproduced as object methods. Their numerical results are available in the
typed Fortran result.

## Share-level simulation

| Upstream numerical component | Fortran implementation | Status |
|---|---|---|
| daily optimization | `simulate_day` | Implemented |
| multi-day lifecycle | `simulate_portfolio` | Implemented |
| corporate-action share adjustment | `adjustment_ratio` input | Implemented |
| internal transfer/market-order split | `allocate_market_and_transfer_orders` | Implemented |
| volume-limited fills | `fill_rate_pct_volume` | Implemented |
| dividends and distributions | `simulate_day` | Implemented |
| transaction costs | `transaction_cost_pct` | Implemented |
| financing costs and financing days | `financing_cost_pct`, `financing_days` | Implemented |
| force trimming | `force_trim_factor` | Implemented |
| force exit of noninvestable positions | `force_exit_non_investable` | Implemented |
| delisting liquidation and return | `delisting`, `delisting_return` | Implemented |
| gross/net P&L and market values | `day_result`, `simulation_result` | Implemented |
| turnover summaries | `simulation_result%turnover` | Implemented |

The upstream orchestration around dates, files, callbacks, persistent result
objects, and report generation is not compiled.

## Exposures and performance

| Upstream computation | Fortran implementation | Status |
|---|---|---|
| `calculate_exposures` factor values | `calculate_exposures` | Implemented |
| category exposure pivots | `calculate_exposures%category` | Implemented with integer level indexing |
| cumulative/annualized returns | `summarize_performance` | Implemented |
| annualized volatility and Sharpe | `summarize_performance` | Implemented |
| maximum drawdown | `summarize_performance` | Implemented |
| average GMV/NMV/turnover | `summarize_performance` | Implemented |
| holding-period estimate | `summarize_performance` | Implemented |

## Excluded infrastructure

- R6 object lifecycle and mutable public/private methods;
- YAML strategy configuration and validation;
- Arrow/Feather and bundled-data readers;
- dplyr/tidyr/tibble/xts-style data manipulation;
- Shiny callbacks and example application;
- plotting, reports, HTML/Word output, and presentation methods;
- file-backed simulation result storage.

These exclusions are interfaces and presentation layers, not additional
portfolio-optimization formulas.
