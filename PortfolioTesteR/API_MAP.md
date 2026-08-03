# API coverage map

This file maps the main upstream `PortfolioTesteR` computational areas to the
Fortran modules. A **direct** entry preserves the core numerical operation. An
**adapted** entry preserves the purpose but replaces R-specific dispatch or data
structures. An **omitted** entry is retained only in the upstream snapshot.

## Data and utilities

| Upstream area | Fortran API | Coverage |
|---|---|---|
| sample price data | `generate_sample_prices` | Adapted deterministic synthetic data |
| panel returns | `panel_returns_simple`, `panel_returns_log` | Direct |
| forward fill | `forward_fill` | Direct |
| alignment | `align_to_indices` | Adapted integer-index alignment |
| invert/standardize signals | `invert_signal`, `standardize_panel` | Direct |
| joins and calendar conversion | none | Omitted; caller aligns arrays explicitly |

## Technical indicators

| Upstream function/family | Fortran API | Coverage |
|---|---|---|
| `calc_momentum` | `calc_momentum` | Direct |
| `calc_distance` | `calc_distance` | Direct |
| moving average | `calc_moving_average` | Direct simple moving average |
| `calc_rsi` | `calc_rsi` | Direct Wilder smoothing |
| stochastic %D | `calc_stochastic_d` | Direct close-only form |
| `calc_cci` | `calc_cci` | Direct close-only form |
| `calc_stochrsi` | `calc_stochrsi` | Direct |
| rolling correlation | `calc_rolling_correlation` | Direct against a selected column |
| rolling volatility | `calc_rolling_volatility` | Direct std/range/MAD/mean-absolute-return; downside extension |
| Bollinger bands | `calc_bollinger_bands` | Direct |
| ATR | `calc_atr` | Adapted close-only true-range proxy |

Volatility methods are selected with `vol_std`, `vol_range`, `vol_mad`,
`vol_abs_return`, or `vol_downside`.

## Filters and selection

| Upstream area | Fortran API | Coverage |
|---|---|---|
| top/bottom N or rank | `filter_top_n` | Direct |
| above/below threshold | `filter_threshold` | Direct |
| between/outside range | `filter_between` | Direct |
| percentile selection | `filter_by_percentile` | Direct |
| combine filters | `combine_filters` | Direct AND/OR |
| regime masks | `apply_regime` | Direct |
| position limits and selection metadata | selection matrices plus `selection_counts` | Adapted |

## Portfolio weighting and risk allocation

| Upstream area | Fortran API | Coverage |
|---|---|---|
| equal weighting | `weight_equally` | Direct |
| signal weighting | `weight_by_signal` | Direct |
| rank weighting | `weight_by_rank` | Direct linear/exponential |
| inverse volatility | `weight_by_volatility` | Direct |
| combine/switch weights | `combine_weights`, `switch_weights` | Direct |
| exposure/group caps | `cap_exposure` | Direct numeric equivalent |
| turnover cap | `cap_turnover` | Direct |
| HRP | `calculate_hrp_weights`, `rolling_hrp_weights` | Native average-linkage implementation |
| risk parity/ERC | `calculate_erc_weights`, `rolling_risk_parity_weights` | Native iterative implementation |
| maximum diversification | `calculate_max_div_weights` | Native projected iteration |

## Backtesting and performance

| Upstream area | Fortran API | Coverage |
|---|---|---|
| `run_backtest` | `run_backtest` | Adapted typed, matrix-first engine |
| portfolio returns | `portfolio_returns_from_weights` | Direct lagged-weight calculation |
| integer/fractional shares | `integer_shares` option | Direct |
| transaction costs | `cost_bps` option | Direct proportional costs |
| stop loss | `stop_loss` option | Adapted same-frequency close check |
| drawdowns | `calculate_drawdown_series` | Direct |
| performance summary | `perf_metrics`, `analyze_performance` | Direct numerical metrics |
| benchmark analysis | `benchmark_statistics` | Direct |
| recovery time | `calculate_recovery_time` | Direct |
| period/regime analysis | `create_regime_buckets` | Adapted numeric buckets |
| print/summary/plot methods | none | Omitted |

## Cross-sectional analytics

| Upstream area | Fortran API | Coverage |
|---|---|---|
| relative-strength ranks | `calc_relative_strength_rank` | Direct |
| top-minus-bottom spreads | `calc_spread_indicators` | Direct |
| correlation dispersion | `calc_correlation_dispersion` | Direct |
| market breadth | `calc_market_breadth` | Direct |
| within-sector/group ranks | `rank_within_groups` | Direct with integer group IDs |
| sector/group breadth | `group_breadth` | Direct |
| group-relative signals | `group_relative_indicators` | Direct |

## Machine-learning helpers

| Upstream area | Fortran API | Coverage |
|---|---|---|
| panel lags/operations | `panel_lag`, `panel_op`, `add_interaction` | Direct |
| labels | `make_labels` | Direct simple/log/sign forms |
| score transforms | `transform_scores` | Direct rank/z-score style forms |
| top-K/group limits | `select_top_k_scores` | Direct |
| score-to-weight mapping | `weight_from_scores` | Direct equal/rank/linear/softmax forms |
| ensemble scores | `combine_scores` | Direct mean/weighted/trimmed forms |
| pooled linear/ridge models | `fit_linear_model`, `predict_linear_model` | Native implementation |
| rolling fit/predict | `rolling_fit_predict` | Native pooled linear/ridge path |
| IC and score evaluation | `ic_series`, `evaluate_scores`, `bucket_returns`, `coverage_by_date` | Direct |
| random forest/XGBoost/deep learning | none | Omitted external engines |
| sequence-model training | none | Omitted external engines |

## Optimization and walk-forward

| Upstream area | Fortran API | Coverage |
|---|---|---|
| parameter-grid search | `run_param_grid` | Adapted procedure-callback interface |
| Sharpe metric | `metric_sharpe` | Direct |
| built-in strategy builders | `momentum_top_n_strategy`, `rsi_reversion_strategy`, `volatility_adjusted_momentum_strategy` | Native examples/reusable callbacks |
| walk-forward splits | `make_walk_forward_splits` | Direct integer-period equivalent |
| walk-forward optimization | `run_walk_forward` | Adapted procedure-callback interface |
| stitched OOS equity | `stitch_returns` | Direct |
| heatmaps/diagnostic plots | none | Omitted |

## R-specific and external infrastructure omitted

- S3 classes, generic methods, formula/list dispatch, and nonstandard evaluation.
- `data.table`, `zoo`, `xts`, dates, calendars, symbol metadata, and names.
- CSV, SQL, Yahoo/Quantmod, web scraping, and SQLite adapters.
- Graphics, reports, vignettes, progress bars, and interactive diagnostics.
- Parallel orchestration and dynamically sourced R strategy functions.
- `glmnet`, `ranger`, `xgboost`, Keras, and TensorFlow model engines.
- Packaged `.rda` datasets; the upstream files remain in `upstream/`.
