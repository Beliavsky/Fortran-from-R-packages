# Translation coverage

## Fully translated computational functionality

- all 14 exported risk/performance-measure wrappers
- all point-estimate formulas used by those wrappers
- generic vector and matrix-column standard-error dispatch
- iid influence-function standard errors
- correlated and AR-prewhitened influence-function standard errors
- adaptive correlated standard errors
- iid and serial block bootstraps
- periodogram generation and frequency selection
- polynomial frequency-domain design construction and optional standardization
- exponential and Gamma elastic-net periodogram fitting through vendored RPEGLMEN
- robust outlier cleaning and all influence functions through vendored RPEIF
- lag-one return and influence-function correlations

## Replaced by documented numerical equivalents

- R `fft`: self-contained direct DFT
- R `arima(..., order=c(1,0,0))`: intercept AR(1) least squares
- R `boot`: internal deterministic iid resampling
- R `tsboot(sim="fixed")`: circular fixed-length block resampling
- R matrix/`xts` dispatch: ordinary rank-1 and rank-2 Fortran arrays

## Omitted as noncomputational or R-specific

- `xts` and `zoo` index preservation
- S3-style list formatting and `printSE`
- plotting and vignette figures
- R documentation, namespace, and package-loading infrastructure
- `PerformanceAnalytics` example-data integration
