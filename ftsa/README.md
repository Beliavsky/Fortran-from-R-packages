# ftsa — modern Fortran translation

This package translates computational routines from R package **ftsa 6.7** to modern free-form Fortran with FPM. It is intended to live as the top-level `ftsa/` directory in [Fortran-from-R-packages](https://github.com/Beliavsky/Fortran-from-R-packages).

The implementation uses `real(dp)` throughout, with `dp` imported from `rfortran-core` (`r_kinds`). Shared numerical work is reused through sibling FPM dependencies rather than copied into this package. In particular, PCA/eigendecomposition operations use the public `rfortran-linalg` API; no system BLAS/LAPACK link flags are required by this package.

## Build

From a repository checkout that also contains sibling directories `rfortran-core/` and `rfortran-linalg/`:

```text
fpm build
fpm test
fpm run --example basic
```

## Implemented numerical areas

The current checkpoint includes regular-grid classical functional time-series PCA, random-walk score forecasting, recursive FTSM forecasts, functional summary statistics and differencing, the full upstream forecast-error metric family, bootstrap means, functional autocorrelation and long-run covariance, ER/GR component selection, rescaled-range estimation, FAR/FARMA simulation, SIMPLS functional regression, multivariate MAF, principal-score bootstrap data, weighted FTSM tuning, rolling ISFE, partial-curve updating, functional linear regression updating, and forward/backward stopping-time monitoring.

## Translation coverage

Package status: **substantial**.

Coverage: **39 of 39 (100.0%)** exported computational R functions.

The fraction measures computational functions with meaningful Fortran mappings; it does **not** mean complete R-interface compatibility. Plotting, printing/presentation, type predicates, generic-only dispatch, and functions that merely dispatch to another package are excluded from the denominator. Formula/S3 behavior, R object construction, exact R RNG streams, and dependency-specific high-level interfaces may be omitted even for mapped functions.

| R function | Fortran module | Fortran procedure(s) | Mapping status |
|---|---|---|---|
| `extract` | `ftsa_statistics` | `ftsa_extract_time`, `ftsa_extract_x` | substantial |
| `ftsm` | `ftsa_models` | `ftsa_ftsm` | partial |
| `forecast.ftsm` | `ftsa_models` | `ftsa_forecast_ftsm` | partial |
| `residuals.fm` | `ftsa_statistics` | `ftsa_residuals` | complete |
| `diff.fts` | `ftsa_statistics` | `ftsa_difference` | substantial |
| `mean.fts` | `ftsa_statistics` | `ftsa_mean_coordinate` | substantial |
| `median.fts` | `ftsa_statistics` | `ftsa_median_coordinate` | substantial |
| `var.fts` | `ftsa_statistics` | `ftsa_variance_coordinate` | substantial |
| `error` | `ftsa_statistics` | `ftsa_error_metric`, `ftsa_error_all` | complete |
| `fbootstrap` | `ftsa_statistics` | `ftsa_bootstrap_mean` | partial |
| `sd.fts` | `ftsa_statistics` | `ftsa_sd_coordinate` | substantial |
| `quantile.fts` | `ftsa_statistics` | `ftsa_quantiles` | substantial |
| `sd.default` | `ftsa_statistics` | `ftsa_sd_vector`, `ftsa_sd_columns` | substantial |
| `ftsmiterativeforecasts` | `ftsa_models` | `ftsa_ftsm_iterative_forecasts` | partial |
| `centre` | `ftsa_statistics` | `ftsa_centre_mean`, `ftsa_centre_variance` | partial |
| `long_run_covariance_estimation` | `ftsa_covariance` | `ftsa_long_run_covariance` | substantial |
| `facf` | `ftsa_covariance` | `ftsa_facf` | complete |
| `ER_GR` | `ftsa_factor` | `ftsa_er_gr` | substantial |
| `stop_time_sim_data` | `ftsa_simulation` | `ftsa_stop_time_sim_data` | substantial |
| `sim_FARMA` | `ftsa_simulation` | `ftsa_sim_farma` | substantial |
| `H_estimate_RS` | `ftsa_factor` | `ftsa_h_estimate_rs` | substantial |
| `isfe.fts` | `ftsa_extended` | `ftsa_isfe` | partial |
| `fplsr` | `ftsa_extended` | `ftsa_fplsr` | substantial |
| `dynupdate` | `ftsa_extended` | `ftsa_dynupdate` | partial |
| `ftsmweightselect` | `ftsa_extended` | `ftsa_ftsm_weight_select` | partial |
| `forecastfplsr` | `ftsa_extended` | `ftsa_forecast_fplsr` | substantial |
| `pcscorebootstrapdata` | `ftsa_extended` | `ftsa_pcscore_bootstrap_data` | partial |
| `MAF_multivariate` | `ftsa_extended` | `ftsa_maf_multivariate` | substantial |
| `dynamic_FLR` | `ftsa_extended` | `ftsa_dynamic_flr` | partial |
| `stop_time_detect` | `ftsa_extended` | `ftsa_stop_time_detect` | partial |
| `MFDM` | `ftsa_remaining` | `ftsa_mfdm` | partial |
| `farforecast` | `ftsa_remaining` | `ftsa_farforecast` | partial |
| `T_stationary` | `ftsa_remaining` | `ftsa_t_stationary` | partial |
| `skew_t_fun` | `ftsa_remaining` | `ftsa_skew_t_fun` | partial |
| `Horta_Ziegelmann_FPCA` | `ftsa_remaining` | `ftsa_horta_ziegelmann_fpca` | partial |
| `CoDa_FPCA` | `ftsa_remaining` | `ftsa_coda_fpca` | partial |
| `LQDT_FPCA` | `ftsa_remaining` | `ftsa_lqdt_fpca` | partial |
| `CoDa_BayesNW` | `ftsa_remaining` | `ftsa_coda_bayes_nw` | partial |
| `GAEVforecast` | `ftsa_gaev` | `ftsa_gaev_forecast` | partial |

All 39 functions in the audited computational denominator have a meaningful Fortran mapping. Many high-level mappings remain intentionally partial where the R implementation relies on package-specific fitting engines, object systems, plotting, or interval/bootstrap interfaces.

See `docs/API_COVERAGE.md` for scope notes and `docs/VALIDATION.md` for the validation record.
