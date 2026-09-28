# API coverage

The coverage basis is the set of distinct **exported computational R functions** in upstream `NAMESPACE`.

Excluded from the denominator are plotting/presentation exports (`plot.ftsm`, `plot.fm`, `plot.fmres`, `plot.ftsf`, `plotfplsr`, `summary.fm`), the type predicate `is.fts`, generic-only exports (`sd`, `var`, `quantile`), and `var.default`, which delegates directly to `stats::var` rather than implementing an ftsa-specific computational algorithm.

The resulting denominator is 39 functions. The current checkpoint maps all 39, for 100.0% mapped-function coverage.

Important compatibility boundaries:

- `ftsm` is the classical regular-grid PCA core. Robust PCA variants, spline interpolation/smoothing, and missing-value reconstruction are not yet translated.
- `forecast.ftsm` currently implements random-walk score forecasts with or without drift. Forecast-package ETS/ARIMA/AR/structural/ARFIMA models and interval object construction remain untranslated.
- `fbootstrap` currently implements nonparametric curve resampling around the mean functional; smoothing and other functional centers remain untranslated.
- Simulation and bootstrap routines preserve the mathematical resampling/DGP structure but intentionally do not claim R RNG-stream identity.
- `fplsr`/`forecastfplsr` currently use the upstream SIMPLS path; NIPALS and prediction-interval branches are not mapped.
- `dynupdate` maps the numerical `ts`, OLS, penalized PLS and ridge updates with random-walk score forecasting; block, ARIMA/ETS and interval branches remain outside the mapping.
- `dynamic_FLR` maps the classical point-forecast eigenspace regression; robust decomposition and maximum-entropy bootstrap intervals are omitted.
- `stop_time_detect` maps the forward/backward error construction but uses a local minimum-SSE single mean-shift selector rather than the external `strucchange` and `ecp` packages.

- `farforecast` maps the functional-PCA/FPE/VAR point-forecast algorithm with a lightly ridge-stabilized score-VAR regression; exact `vars`/`arima` internals and interval objects are omitted.
- `T_stationary` maps the long-run covariance statistic, pivotal/nonpivotal variants, Fourier projection, and Monte Carlo calibration; analytic Fourier functions replace `fda` objects.
- `Horta_Ziegelmann_FPCA`, `CoDa_FPCA`, `LQDT_FPCA`, and `CoDa_BayesNW` map their density-transform/eigenspace/regression/forecast kernels for already discretized densities; raw-sample KDE and several package-specific bootstrap/model-selection branches are omitted.
- `skew_t_fun` uses the same Fernandez-Steel standardized skew-Student density formulas as the upstream fGarch dependency, with log-density likelihood evaluation plus a stabilized local VAR/AIC forecast layer.
- `MFDM` maps the common and sex-specific two-level SVD decomposition and deterministic score forecasts; the central JAGS posterior-sampling and posterior-predictive uncertainty layer remains outside the Fortran implementation.
- `GAEVforecast` maps smooth GEV fitting, cross-validated basis dimensions, ridge-stabilized trend/VAR coefficient forecasting, mode reconstruction, and density evaluation; its Fortran-native penalized open-uniform B-spline fit is not an exact reproduction of `evgam` cubic-regression-spline smoothing or smoothing-parameter selection.
- `facf` returns zero at lag zero because the Fortran numeric result array has no R-style NA sentinel contract; all nonzero-lag computations follow the upstream Frobenius-norm definition.
- R plotting, formulas, S3 classes, `fts`/`ts` metadata and name handling are outside the numerical API.

The authoritative per-function mapping metadata is in `fpm.toml`.
