# Translation coverage

## Computational exports

| R export | Fortran procedure |
|---|---|
| `ActivePremium` | `active_premium` |
| `AdjustedSharpeRatio` | `adjusted_sharpe_ratio` |
| `BernardoLedoitRatio` | `bernardo_ledoit_ratio` |
| `BurkeRatio` | `burke_ratio` |
| `DRatio` | `d_ratio` |
| `KellyRatio` | `kelly_ratio` |
| `MartinRatio` | `martin_ratio` |
| `SkewnessKurtosisRatio` | `skewness_kurtosis_ratio` |
| `PainIndex` | `pain_index` |
| `MeanAbsoluteDeviation` | `mean_absolute_deviation` |
| `CalmarRatio` | `calmar_ratio` |
| `SterlingRatio` | `sterling_ratio` |
| `AppraisalRatio` | `appraisal_ratio` |
| `TrackingError` | `tracking_error` |
| `InformationRatio` | `information_ratio` |
| `TreynorRatio` | `treynor_ratio` |
| `DownsideDeviation` | `downside_deviation` |
| `OmegaSharpeRatio` | `omega_sharpe_ratio` |
| `SortinoRatio` | `sortino_ratio` |
| `ProspectRatio` | `prospect_ratio` |
| `VolatilitySkewness` | `volatility_skewness` |
| `M2Sortino` | `m2_sortino` |
| `SharpeRatio` | `sharpe_ratio` |
| `SharpeRatio.annualized` | `sharpe_ratio_annualized` |
| `PainRatio` | `pain_ratio` |
| `table.AnnualizedReturns` | `table_annualized_returns` |
| `Return.annualized` | `return_annualized` |
| `CAPM.jensenAlpha` | `capm_jensen_alpha` |
| `UlcerIndex` | `ulcer_index` |
| `DrawdownPeak` | `drawdown_peak` |
| `maxDrawdown` | `max_drawdown` |
| `durbinH` | `durbin_h` |

The internal R helpers for skewness, kurtosis, upside/downside risk, historical
VaR/ES, and drawdown paths are also exposed as typed Fortran procedures.

## Not compiled

- `getEER`
- `getFed`
- `getFrench.Factors`
- `getFrench.Portfolios`
- plotting and `xts`/data-frame presentation infrastructure
- bundled R data objects
