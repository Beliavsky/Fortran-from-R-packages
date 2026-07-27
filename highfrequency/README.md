# highfrequency-fortran

Modern Fortran 2018/FPM computational port of the numerical core of the R
package `highfrequency` 1.0.2.

The upstream package provides tools for high-frequency trade and quote data,
realized volatility and covariance, jump detection, lead-lag estimation,
liquidity measurement, spot volatility and drift, and HAR/HEAVY forecasting.
This port focuses on self-contained numerical algorithms operating on ordinary
Fortran arrays. It does not attempt to reproduce R's `xts`, `data.table`,
formula, plotting, network-download, or serialization infrastructure.

## License

The upstream package declares `GPL (>= 2)`. This derivative project is
distributed under `GPL-2.0-or-later`. The complete GPL version 2 text is in
`LICENSE`. The original package and its attribution are retained under
`original/highfrequency-1.0.2`.

## Build

```text
fpm build
fpm test
fpm run highfrequency_demo
fpm run --example realized_measures
fpm run --example microstructure
```

No external numerical library is required.

## Data conventions

- Observation times are integer ticks chosen by the caller, such as seconds or
  milliseconds from an origin.
- Return matrices have observations in rows and assets in columns.
- Price vectors must be positive for logarithmic-return routines.
- Missing values are not encoded as R-style `NA`; clean or mask observations
  before calling a numerical routine.
- `make_returns` returns `n-1` returns from `n` prices. Upstream R
  `makeReturns` retains length `n` and inserts a leading zero. The Fortran
  convention avoids a synthetic observation.

## Main modules

The umbrella module is:

```fortran
use highfrequency
```

Important interfaces include:

```fortran
returns = make_returns(prices)
rv = rrvar(returns)
cov = rcov(return_matrix)
bpv = rbpvar(returns)
cov = rkernelcov(return_matrix, 10, "parzen", force_psd=.true.)
hy = rhy_cov(times_x, prices_x, times_y, prices_y)

har = fit_har(daily_realized_variance, [1, 5, 22])
forecast = har_forecast(har, daily_realized_variance, [1, 5, 22])

heavy = fit_heavy(daily_returns, daily_realized_measure)
jump = bns_jump_test(intraday_returns)

call liquidity_measures(price, size, bid, offer, bid_size, offer_size, &
  window, liquidity)
```

Both descriptive names and familiar package-style aliases are provided for
many realized measures, including `rrvar`, `rcov`, `rbpvar`, `rbpcov`,
`rskew`, `rkurt`, `rsvar`, `rquar`, `rtpquar`, `rqpvar`, `rminrvar`,
`rmedrvar`, `rminrquar`, `rmedrquar`, `rmpvar`, `rthresholdcov`,
`rkernelcov`, `rtsvar`, `rtscov`, `rhy_cov`, and `rbeta`.

## Implemented areas

- Log and arithmetic return construction
- Fixed-time, business-time, and OHLCV aggregation primitives
- Previous-tick and pairwise refresh-time synchronization
- Trade/quote matching and cleanup masks
- Trade-direction and liquidity diagnostics
- Realized variance, covariance, skewness, kurtosis, semivariance, and
  semicovariance
- Bipower, tripower, quadpower, minimum, median, multipower, threshold,
  pre-averaged, average, kernel, modulated, and two-scale estimators
- Hayashi-Yoshida covariance and realized beta
- ReMeDI covariance sequence and automatic `kn` selection
- BNS, ABD, and Ait-Sahalia-Jacod jump tests
- Integrated-variance confidence intervals
- Lead-lag contrast estimation
- Kernel spot volatility, spot drift, and drift-burst statistics
- Basic HAR regression and forecasts
- Two-equation HEAVY quasi-maximum-likelihood fitting and forecasts
- Positive-semidefinite projection, covariance-to-correlation conversion,
  OLS, bounded Nelder-Mead, and supporting statistics

See `COVERAGE.md` for routine-level scope and exclusions.

## Important scope limits

The following upstream facilities are retained for reference but are not
compiled into the Fortran library:

- `xts`, `zoo`, `data.table`, R formulas, S3 methods, and multi-day dispatch
- Plotting and graphics
- Alpha Vantage downloads
- RDS/RDA data loading and package data objects
- Exchange-code and sales-condition tables tied to R data sets
- Exact reproductions of robustbase MCD, Rsolnp, sandwich, and numDeriv
  adapters
- Several specialized multivariate estimators and rank/bootstrap jump tests
  listed in `COVERAGE.md`

The translated numerical routines are intended as a tested foundation that can
be extended without requiring an R runtime.

## Repository layout

```text
src/        library modules
test/       FPM tests
app/        demonstration program
example/    focused examples
scripts/    direct compiler validation scripts
original/   unmodified upstream package
provenance/ supplied archive and checksum manifests
```
