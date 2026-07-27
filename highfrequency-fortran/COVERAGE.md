# Computational coverage

Source package: `highfrequency` 1.0.2  
Upstream license: GPL version 2 or later

This project translates the self-contained numerical core into array-based
Fortran. R container dispatch, presentation, network access, and data-set
lookups are outside the compiled scope.

## Fully implemented numerical families

| Upstream area | Fortran interfaces |
|---|---|
| `makeReturns` | `make_returns` |
| `aggregatePrice`, `aggregateTS` numerical core | `aggregate_last`, `aggregate_sum` |
| `businessTimeAggregation` numerical grouping | `business_time_groups` |
| `makeOHLCV` | `make_ohlcv` |
| `mergeTradesSameTimestamp`, `mergeQuotesSameTimestamp` numerical core | `merge_same_timestamp` |
| `matchTradesQuotes` previous-tick mode | `match_trades_quotes` |
| `refreshTime` pairwise core | `refresh_time_pair` |
| `noZeroPrices`, `noZeroQuotes` | `no_zero_prices_mask`, `no_zero_quotes_mask` |
| `rmNegativeSpread`, `rmLargeSpread` | `nonnegative_spread_mask`, `maximum_spread_mask` |
| rolling price outlier cleanup | `price_outlier_mask` |
| `spreadPrices` | `spread_prices` |
| `getTradeDirection` | `trade_direction` |
| `getLiquidityMeasures` principal measures | `liquidity_measures` |
| `RV`, `rRVar` | `realized_variance`, `rrvar` |
| `rCov` | `realized_covariance`, `rcov` |
| `rSkew`, `rKurt` | `realized_skewness`, `realized_kurtosis`, `rskew`, `rkurt` |
| `rSV`, `rSVar` | `realized_semivariance`, `rsvar` |
| `rSemiCov` | `realized_semicovariance`, `rsemicov` |
| `rBPCov` | `bipower_variation`, `bipower_covariance`, `rbpvar`, `rbpcov` |
| `rQuar` | `realized_quarticity`, `rquar` |
| `rTPQuar` | `tripower_quarticity`, `rtpquar` |
| `rQPVar` | `quadpower_variation`, `rqpvar` |
| `rMinRVar`, `rMedRVar` | `minimum_realized_variance`, `median_realized_variance` |
| `rMinRQuar`, `rMedRQuar` | `minimum_realized_quarticity`, `median_realized_quarticity` |
| `rMPVar` | `multipower_variation`, `rmpvar` |
| `rThresholdCov` | `threshold_covariance`, `rthresholdcov` |
| `rKernelCov` principal kernels | `realized_kernel_covariance`, `rkernelcov` |
| pre-averaging estimators | `preaverage_returns`, `preaveraged_covariance` |
| `rAVGCov` synchronous core | `average_realized_covariance` |
| `rMRC`, `rMRCov` numerical core | `modulated_realized_covariance` |
| `rTSCov` synchronous pairwise core | `two_scale_variance`, `two_scale_covariance`, `rtsvar`, `rtscov` |
| `rHYCov` | `hayashi_yoshida_covariance`, `rhy_cov` |
| `rBeta` | `realized_beta`, `rbeta` |
| microstructure-noise estimate | `noise_variance` |
| `ReMeDI` | `remedi` |
| `knChooseReMeDI` | `choose_remedi_kn` |
| `BNSjumpTest` | `bns_jump_test` |
| `ABDJumptest` | `abd_jump_test` |
| `AJjumpTest` | `aj_jump_test` |
| `IVinference` | `iv_inference` |
| `leadLag` | `lead_lag` |
| `spotVol` numerical kernel estimator | `spot_volatility` |
| `spotDrift` mean and median estimators | `spot_drift` |
| `driftBursts` core statistic | `drift_burst_statistic` |
| basic `HARmodel` | `fit_har`, `har_forecast` |
| `HEAVYmodel` | `fit_heavy`, `heavy_recursion`, `heavy_forecast` |
| `makePsd` | `make_psd` |

## Partial compatibility

| Upstream routine/family | Status |
|---|---|
| `aggregateTrades`, `aggregateQuotes` | Aggregation primitives are implemented; R column-name and per-symbol dispatch are not. |
| `tradesCleanup`, `quotesCleanup`, `tradesCleanupUsingQuotes` | Numerical masks and matching are implemented; exchange calendars and R condition-code tables are not. |
| `matchTradesQuotes` | Previous-tick matching is implemented. The optional backwards-forwards matching report and graphics are excluded. |
| `getLiquidityMeasures` | Core spreads, impact, imbalance, size, and proportional measures are implemented. R table columns and lagged return presentation are not reproduced as an S3 object. |
| `HARmodel` | Standard HAR with arbitrary averaging periods, optional log transform, and external regressors is implemented. Specialized HARJ/HARQ/CHAR formula construction and R inference summaries are not. |
| `HEAVYmodel` | Both recursions, constrained QML estimation, and forecasting are implemented. Sandwich standard errors and R plotting/summary methods are not. |
| `spotVol`, `spotDrift`, `driftBursts` | Core kernel estimators are implemented. R's complete collection of pre-averaging, periodicity, and multi-day adapters is not. |
| `rKernelCov` | Bartlett, Parzen, Tukey-Hanning, and flat-top kernels are implemented. The complete upstream kernel catalogue and automatic bandwidth selectors are not. |
| `rTSCov`, `rRTSCov` | Synchronous pairwise two-scale estimation is implemented. Full multivariate refresh-time assembly and robust tuning are not. |
| `rMRCov` | A self-contained pre-averaged/modulated estimator is implemented, not an exact RcppArmadillo reproduction of every tuning path. |

## Not compiled

These routines depend substantially on R infrastructure, bundled data,
external packages, or specialized algorithms not yet translated:

- `getAlphaVantageData`
- `plotTQData` and all `plot.*`, `print.*`, and `summary.*` methods
- exchange/sales-condition selection helpers tied to R tables
- `rOWCov` exact robustbase MCD implementation
- `rCholCov` permutation-selection framework
- `rBACov`
- `rRTSCov` full robust multivariate implementation
- `ReMeDIAsymptoticVariance`
- `JOjumpTest`
- `intradayJumpTest` complete Lee-Mykland object workflow
- `rankJumpTest` bootstrap/SVD workflow
- `getCriticalValues` and bundled critical-value data
- exact drift-burst parallel/Rcpp implementation
- R-specific `makeRMFormat`, `gatherPrices`, and multi-day `xts` dispatch

The original implementations remain under `original/highfrequency-1.0.2`.
