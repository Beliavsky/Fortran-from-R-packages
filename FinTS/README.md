# FinTS-fortran

A modern Fortran 2018/FPM translation of the computational routines in the R package **FinTS 0.4-9**, the companion package to Ruey S. Tsay's *Analysis of Financial Time Series*.

## Scope

The port includes the numerical content of the exported functions:

- `Acf`: univariate autocorrelation, autocovariance, and partial autocorrelation
- multivariate cross-autocorrelation through the additional `cross_acf` routine
- `ArchTest`: Engle ARCH LM test
- `AutocorTest`: Ljung-Box, Box-Pierce, and rank Ljung-Box tests
- `ARIMA`: nonseasonal and seasonal ARIMA fitting with optional regressors and residual Box tests
- `FinTS.stats`: count, mean, standard deviation, skewness, excess kurtosis, minimum, and maximum
- `apca`: asymptotic principal component analysis
- `plotArmaTrueacf`: computational core only, including AR roots, stationarity, theoretical ACF/PACF, damping, and periodicity
- `findConjugates`: complex-conjugate-pair detection
- `compoundInterest` and `simple2logReturns`
- `as.yearmon2`: numeric `yyyy.mm` and integer `yyyymm` conversion with duplicate-month detection

The following R runtime or data-management facilities are not included:

- plotting methods and graphics
- `package.dir`
- `runscript`
- `url2data`
- the general `read.yearmon` wrapper around `read.table` and `zoo`
- bundled data sets and chapter scripts
- R S3 classes, `zoo` indices, and data-frame formatting

Original metadata, documentation, and R source are retained under `original/` for attribution and provenance.

## Build with FPM

```console
fpm build
fpm test
fpm run
fpm run --example time_series_diagnostics
fpm run --example arima_example
fpm run --example apca_example
```

No external Fortran library is required.

A direct GNU Fortran validation script is also included:

```console
./run_gfortran_tests.sh
```

On Windows:

```console
run_gfortran_tests.bat
```

## Minimal example

```fortran
program example
   use fints
   implicit none

   real(dp) :: x(200)
   type(arima_result) :: fit
   integer :: i

   x(1) = 2.0_dp
   do i = 2, size(x)
      x(i) = 2.0_dp + 0.6_dp * (x(i - 1) - 2.0_dp) + &
         0.1_dp * sin(0.17_dp * real(i * i, dp))
   end do

   call ARIMA(x, [1, 0, 0], fit)
   print '(a,f10.5)', 'mean = ', fit%intercept
   print '(a,f10.5)', 'AR(1) = ', fit%ar(1)
   print '(a,f10.5)', 'Box p-value = ', fit%box_test%p_value
end program example
```

## API conventions

R lists and test objects are represented by derived types:

- `acf_result`
- `cross_acf_result`
- `test_result`
- `summary_result`
- `apca_result`
- `arma_acf_result`
- `yearmon_result`
- `arima_result`

Status codes are available as `fints_ok`, `fints_invalid_input`, `fints_singular`, `fints_nonstationary`, `fints_iteration_limit`, `fints_no_data`, `fints_numerical_failure`, and `fints_io_error`.

## ARIMA implementation

The R function delegates fitting to `stats::arima`, which includes exact likelihood and Kalman filtering. This dependency-free port uses conditional Gaussian likelihood and a Nelder-Mead optimizer. It supports:

- ARIMA `(p,d,q)` orders
- seasonal `(P,D,Q)` orders and an explicit seasonal period
- optional regressors
- optional mean for undifferenced models
- reflection-coefficient transformations for stationary AR and invertible MA parameters
- `CSS`, `ML`, and `CSS-ML` method labels
- conditional residuals, Gaussian log likelihood, AIC, and residual Box tests

All three method labels currently use the same conditional Gaussian objective. Numerical estimates can therefore differ from R, especially for short samples, differenced models, missing data, or models near a stationarity boundary. See `PORTING.md`.

## License

The original package declares `GPL (>= 2)`. This translation is distributed under **GPL-2.0-or-later**. See `LICENSE` and `NOTICE.md`.
