# Upstream API map

| R export | Fortran computational API | Notes |
|---|---|---|
| `DHSimulate` | `dh_simulate`, `dh_condition`, `dhsimulate` | Correct zero-mean Davies-Harte simulation by default; source endpoint behavior optional. |
| `DLAcfToAR` | `dl_acf_to_ar`, `dlacftoar` | Returns `type(dl_ar_result)`. |
| `DLLoglikelihood` | `dl_loglikelihood`, `dlloglikelihood` | Concentrated exact likelihood. |
| `DLResiduals` | `dl_residuals`, `dlresiduals` | Raw or standardized one-step errors. |
| `DLSimulate` | `dl_simulate`, `dlsimulate` | Optional explicit innovation vector replaces R's `rand.gen`. |
| `exactLoglikelihood` | `exact_loglikelihood`, `exactloglikelihood` | Returns exact likelihood and variance estimate. |
| `PredictionVariance` | `prediction_variance`, `predictionvariance` | Durbin-Levinson approximation or exact Toeplitz calculation. |
| `innovationVariance` | `innovation_variance`, `innovationvariance` | AR/AIC and Kolmogoroff periodogram methods. |
| `SimGLP` | `sim_glp`, `simglp` | Direct finite MA convolution. |
| `ToeplitzInverseUpdate` | `toeplitz_inverse_update`, `toeplitzinverseupdate` | Bordering update. |
| `TrenchForecast` | `trench_forecast`, `trenchforecast` | Returns `type(forecast_result)`. |
| `TrenchInverse` | `trench_inverse`, `trenchinverse` | Trench recursion with checked SPD fallback. |
| `TrenchLoglikelihood` | `trench_loglikelihood`, `trenchloglikelihood` | Concentrated likelihood. |
| `TrenchMean` | `trench_mean`, `trenchmean` | Exact GLS/BLUE mean. |
| `is.toeplitz` | `is_toeplitz`, `istoeplitz` | Tolerance-aware symmetric Toeplitz test. |
| `tacvfARMA` | `tacvf_arma`, `tacvfarma` | ARMA autocovariances and stationarity validation. |

Supporting public routines include `toeplitz_matrix`, `ar_to_ma`,
`ar_is_stationary`, `set_ltsa_seed`, `ltsa_normal`, and status constants.
