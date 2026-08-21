# TSA 1.3.1 API map

This map follows the exports in the original TSA `NAMESPACE`.

| TSA/R routine | Fortran routine | Status / notes |
|---|---|---|
| `ARMAspec` | `arma_spectrum` | Implemented computational spectrum, including seasonal AR/MA factors. |
| `BoxCox.ar` | `boxcox_ar` | Implemented AR likelihood profile, MLE grid choice, and 95% profile interval. |
| `Keenan.test` | `keenan_test` | Implemented. |
| `LB.test` | `lb_test` | Implemented Ljung-Box and Box-Pierce variants. |
| `McLeod.Li.test` | `mcleod_li_test` | Implemented. |
| `Tsay.test` | `tsay_test` | Implemented quadratic-regression nonlinearity test. |
| `acf` | `autocorrelation`, `autocovariance`, `partial_autocorrelation`, `cross_correlation` | Computational ACF/PACF/CCF implemented; R plotting/ts metadata omitted. |
| `arima` | `arima_fit` | CSS/CSS-ML plus diffuse state-space ML implemented, including missing-observation updates, Gardner/AS154 stationary initialization, fixed/init NaN masks, AR/SAR stationarity transforms, MA/SMA invertibility normalization, SVD-conditioned all-free xreg, and R-style Hessian covariance. |
| `arima.boot` | `arima_bootstrap`, `arima_bootstrap_sample` | Regular ARIMA bootstrap translated including TSA conditional prefix behavior, Gaussian/residual resampling, and the source MA-convolution convention. |
| `arimax` | `arimax_fit`, `transfer_filter`, `io_regressor` | Ordinary regressors, rational `xtransf` transfer blocks, `io` pulse interventions, fixed/init masks, ARIMA/seasonal errors, SVD regression conditioning, and covariance back-transforms are implemented jointly. |
| `armasubsets` | `armasubsets_fit` | Implemented using the supplied `leaps` Fortran translation. |
| `detectAO` | `detect_ao` | Chang-Chen-Tiao additive-outlier statistic implemented. |
| `detectIO` | `detect_io` | Innovational-outlier statistic and robust scale implemented. |
| `eacf` | `eacf` | Chan extended ACF recursion and coded significance table data implemented. |
| `gBox` | `gbox_test` | Implemented against supplied `tseries` `garch_result`. |
| `garch.sim` | `garch_sim` | Implemented for Gaussian innovations. R's arbitrary RNG callback is not part of the Fortran API. |
| `harmonic` | `harmonic_matrix` | Implemented. |
| `kurtosis` | `kurtosis` | Implemented using TSA's population-moment convention; returns excess kurtosis. |
| `lagplot` | — | Omitted: plotting routine; its locfit/mgcv dependencies are therefore not compiled. |
| `periodogram` | `periodogram` | Raw one-sided DFT periodogram with TSA doubling convention implemented. |
| `prewhiten` | `prewhiten_filter` | Computational prewhitening filter implemented. |
| `qar.sim` | `qar_sim` | Implemented. |
| `runs` | `runs` | Exact run count/distribution calculation implemented. |
| `season` | `season_index` | Implemented numeric seasonal index generation. |
| `skewness` | `skewness` | Implemented using TSA's population-SD convention. |
| `spec` | `spec_pgram`, `spec_ar`, `periodogram`, `arma_spectrum` | Computational `spec.pgram` includes univariate/multivariate spectra, coherence/phase, scalar/per-series tapering, compact/full custom kernels, and named Daniell/modified-Daniell/Fejer/Dirichlet kernels. `spec_ar` supports Yule-Walker, Burg/Burg2, OLS, and ML with fixed order or AIC selection. Plotting is omitted. |
| `tar` | `tar_fit`, `tar_fit_multi` | Univariate and matrix/multiseries SETAR implemented with MAIC/CLS threshold search, optional order selection, missing-window omission, transforms, and per-series diagnostics. |
| `tar.sim` | `tar_sim` | Implemented. |
| `tar.skeleton` | `tar_skeleton` | Implemented deterministic skeleton/cycle search. |
| `tlrt` | `tlrt_test`, `tlrt_p_value` | Implemented using Chan (1990) approximation encoded by TSA. |
| `zlag` | `lag_vector` | Implemented with caller-selected padding value (zero by default). |
| `predict.TAR` | `tar_predict` | Simulation median and 95% interval implemented. |
| `fitted.*`, `rstandard.*` | fields in `arimax_result` / `tar_result` | Numerical results are stored directly; S3 methods are unnecessary. |
| `tsdiag.*` | `lb_test`, `gbox_test`, residual result fields | Computational diagnostics available; graphical diagnostics omitted. |
| `plot.*`, `summary.armasubsets` | — | Presentation-only code omitted. |
