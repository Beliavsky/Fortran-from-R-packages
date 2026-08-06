# Translation coverage

## Complete computational coverage

- general estimator dispatcher
- nuisance-parameter construction
- shape-grid construction
- mean influence function
- robust mean influence function
- SD and semi-SD influence functions
- VaR and ES influence functions
- LPM influence functions, orders 1 and 2
- Sharpe, Sortino, and downside Sharpe influence functions
- ES-ratio and VaR-ratio influence functions
- Omega and Rachev-ratio influence functions
- robust winsorizing filter
- AR(p) prewhitening
- normal PDF, CDF, and inverse CDF
- R type-7 quantiles
- direct Gaussian KDE at a point

## Deliberately omitted

- all graphics and plotting
- `xts`/`zoo` indexes and date restoration
- R argument-class checks that do not apply to typed Fortran arrays
- R namespace and S3 infrastructure
- vignette rendering and R documentation machinery

## Equivalent rather than exact

- direct KDE evaluation replaces the FFT/interpolation implementation in R
- conditional OLS AR(p) replaces `stats::arima`
- robust location iteration uses the translated RobStatTM numerical kernels but
  not the R object wrapper
