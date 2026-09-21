# Time Series Analysis

This is an independent, filtered adaptation of the
[CRAN Task View: Time Series Analysis](https://CRAN.R-project.org/view=TimeSeries),
maintained by Rob J Hyndman and Rebecca Killick, version 2026-09-21. The
[source task view](https://github.com/cran-task-views/TimeSeries/blob/main/TimeSeries.md)
provides the broader annotated guide to R packages.

This page includes 39 translated packages from this repository. Its summaries
are original descriptions of the high-level computational capabilities present
in the Fortran translations, rather than copies of the CRAN annotations.
Plotting, interactive displays, R time-series classes, and other untranslated
R infrastructure are not included. These translations are experimental; see
each package's documentation for validation status and precise API coverage.

## General analysis and diagnostics

- [`roll`](../roll/) computes rolling and expanding summaries, quantiles,
  covariance, correlation, and weighted linear regression.
- [`tseries`](../tseries/) provides ARMA and GARCH fitting, stationarity and
  cointegration tests, nonlinear-dependence diagnostics, resampling methods,
  surrogate generation, and selected financial time-series calculations.
- [`tsa`](../tsa/) provides time-series summaries, correlation and spectral
  analysis, ARMA and seasonal ARIMA modeling, threshold autoregression,
  simulation, forecasting, and model diagnostics.
- [`FinTS`](../FinTS/) provides correlation and partial-correlation analysis,
  ARCH and serial-correlation tests, ARIMA modeling, distribution summaries,
  asymptotic principal components, and financial-return transformations.
- [`timsac`](../timsac/) provides correlation, spectral estimation,
  multivariate time-series analysis, matrix filtering, and related numerical
  routines from the TIMSAC collection.

## Times and dates

- [`lubridate`](../lubridate/) provides typed date and fixed-offset date-time
  parsing, calendar arithmetic, intervals, periods, durations, accessors, and
  rounding operations useful in time-series workflows.

## Forecasting and univariate modeling

- [`imputeTS`](../imputeTS/) imputes missing time-series values using
  interpolation, smoothing, seasonal methods, and model-based estimation.
- [`fracdiff`](../fracdiff/) estimates, simulates, differences, and forecasts
  fractionally integrated ARMA processes and provides long-memory estimators
  and model diagnostics.
- [`arfima`](../arfima/) fits, simulates, and forecasts ARFIMA and related
  long-memory models with seasonal differencing, regressors, and dynamic
  transfer functions.
- [`ltsa`](../ltsa/) provides computational methods for stationary and
  long-memory time-series analysis, including autocovariance recursions and
  Toeplitz-system calculations.
- [`greybox`](../greybox/) provides regression-based forecasting, dynamic and
  ARIMA-error models, model selection and combination, rolling evaluation, and
  intermittent-demand methods.
- [`smoots`](../smoots/) estimates smooth trends and derivatives, selects
  smoothing bandwidths, models ARMA residual structure, and constructs
  analytical and bootstrap forecasts.
- [`forecast`](../forecast/) provides forecasting methods for
  exponential-smoothing, ARIMA, BATS, TBATS, and related time-series models.
- [`KFAS`](../KFAS/) provides Kalman filtering and smoothing for Gaussian and
  exponential-family state-space models, including diffuse initialization,
  non-Gaussian approximation, and importance sampling.
- [`dlm`](../dlm/) provides dynamic-linear-model construction, Kalman
  filtering and smoothing, likelihood evaluation, forecasting, simulation,
  maximum-likelihood fitting, and selected Bayesian sampling methods.

## Structural changes

- [`ecp`](../ecp/) detects multiple change points in univariate and multivariate
  data using energy statistics, hierarchical segmentation, and kernel methods.
- [`strucchange`](../strucchange/) provides structural-change tests,
  breakpoint estimation, fluctuation processes, confidence intervals, and
  monitoring boundaries for regression and time-series models.
- [`changepoint`](../changepoint/) detects changes in means, variances,
  distributions, and regression relationships using AMOC, PELT, binary
  segmentation, segment-neighborhood, CUSUM, and CROPS methods.

## Conditional variance models

- [`rugarch`](../rugarch/) provides filtering, likelihood estimation,
  simulation, forecasting, diagnostics, and risk calculations for a broad
  collection of univariate GARCH and ARFIMA-GARCH specifications.
- [`tsgarch`](../tsgarch/) provides filtering and estimation for univariate
  GARCH-family volatility models with multiple conditional innovation
  distributions.
- [`fGarch`](../fGarch/) provides GARCH, APARCH, and EGARCH simulation and
  estimation, conditional distributions, volatility forecasts, diagnostics,
  and parametric risk measures.

## Frequency analysis, decomposition, and filtering

- [`kza`](../kza/) provides Kolmogorov-Zurbenko and adaptive filtering,
  local-variance estimation, and spectral calculations.
- [`wavelets`](../wavelets/) computes discrete and maximal-overlap wavelet
  transforms, their inverses, and multiresolution analyses.
- [`waveslim`](../waveslim/) provides discrete, maximal-overlap, packet,
  multidimensional, and dual-tree wavelet transforms together with
  multiresolution analysis, denoising, wavelet dependence measures, and
  long-memory spectral methods.
- [`wavethresh`](../wavethresh/) provides discrete and stationary wavelet
  transforms, wavelet packets, thresholding, shrinkage, density estimation,
  spectral analysis, and multidimensional wavelet methods.
- [`signal`](../signal/) provides digital-filter design and application,
  frequency and impulse responses, interpolation and resampling,
  Savitzky-Golay filtering, windows, chirps, and spectrogram calculations.
- [`Rssa`](../Rssa/) provides singular spectrum analysis, reconstruction,
  forecasting, gap filling, grouping, and multivariate and oblique variants.
- [`FKF`](../FKF/) provides fast multivariate linear-Gaussian Kalman filtering
  and smoothing, including support for missing observations.

## Stationarity, cointegration, and nonlinear dynamics

- [`tsdyn`](../tsdyn/) fits and simulates linear and nonlinear autoregressive,
  threshold, smooth-transition, VAR, VECM, TVAR, and TVECM models and provides
  forecasting, impulse-response, bootstrap, and specification-test methods.
- [`tserieschaos`](../tserieschaos/) provides delay embedding, mutual
  information, false-nearest-neighbor analysis, correlation integrals,
  recurrence analysis, Lyapunov estimation, chaotic systems, and numerical
  trajectory integration.
- [`fnonlinear`](../fnonlinear/) provides chaotic maps and systems, nonlinear
  time-series embeddings and dependence measures, recurrence and Lyapunov
  calculations, and tests for nonlinear structure.
- [`urca`](../urca/) performs unit-root, stationarity, cointegration, and
  structural-break analysis.

## Multivariate time-series models

- [`varmapack`](../varmapack/) constructs and simulates VAR, VMA, VARMA, and
  VARMAX models, with sample and theoretical autocovariances and impulse
  responses.
- [`MTS`](../MTS/) provides VAR, VARMA, VMA, VARX, VECM, factor, and
  multivariate volatility models together with forecasting, impulse responses,
  dependence diagnostics, and missing-observation estimation.
- [`tsmarch`](../tsmarch/) provides DCC, asymmetric DCC, copula-GARCH, and
  GO-GARCH workflows with multivariate volatility forecasting, higher-order
  conditional moments, diagnostics, simulation, and risk calculations.
- [`vars`](../vars/) provides VAR estimation and lag selection, restrictions,
  forecasting, impulse responses, forecast-error variance decompositions,
  causality tests, and structural VAR identification.
- [`MARSS`](../MARSS/) provides multivariate autoregressive state-space model
  construction, Kalman filtering and smoothing, likelihood and EM fitting,
  forecasting, simulation, uncertainty calculations, and diagnostics.

## Continuous-time models

- [`sde`](../sde/) provides simulation, transition approximations, parameter
  inference, estimating functions, and nonparametric methods for stochastic
  differential equations and diffusion processes.

## Resampling

- [`boot`](../boot/) provides general nonparametric and weighted bootstrap
  sampling, block bootstrap methods, bootstrap confidence intervals,
  empirical influence calculations, importance weighting, tilting, and
  censored-data resampling primitives.

Some packages appear under only their principal topic even when they provide
methods relevant to several sections. Consult the package README for its full
translated computational scope.
