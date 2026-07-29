# vasicekfit-fortran

A dependency-free modern Fortran/FPM implementation of the numerical algorithms
in the R package `vasicekfit` 0.2.0.

The library fits the extended Vasicek single-factor credit-loss model in which
baseline probability of default is modified by macroeconomic covariates. The
closed-form estimator is an ordinary least-squares regression after transforming
observed loss rates to probit space.

## Features

- Vasicek density, CDF, quantile, and random generation
- Optional macro-factor sensitivities in the loss distribution
- Closed-form probit-OLS estimation of baseline PD, asset correlation, and
  macro-factor sensitivities
- Optional finite-portfolio correction for observed rates containing zero or one
- Optional small-sample variance bias correction
- Fitted values and residuals in probit space
- Conditional mean and conditional quantile prediction
- IID delta-method covariance and Wald confidence intervals
- Self-contained Bartlett/Newey-West HAC covariance with optional lag selection
- Typed results and explicit error messages
- No BLAS, LAPACK, statistics, or optimization-library dependency

## Build with FPM

```text
fpm build
fpm test
fpm run vasicekfit_demo
fpm run --example distribution_functions
fpm run --example fit_and_predict
```

## Basic fitting example

```fortran
use vasicekfit, only : dp, vasicek_fit_result, fit_vasicek

real(dp) :: loss_rate(8), macro(8,1)
type(vasicek_fit_result) :: fit

loss_rate = [0.012_dp, 0.018_dp, 0.021_dp, 0.028_dp, &
             0.035_dp, 0.047_dp, 0.061_dp, 0.072_dp]
macro(:,1) = [-1.4_dp, -1.0_dp, -0.6_dp, -0.2_dp, &
               0.2_dp,  0.6_dp,  1.0_dp,  1.4_dp]

fit = fit_vasicek(loss_rate, macro)
if (.not. fit%ok) error stop fit%message

print *, fit%p
print *, fit%rho
print *, fit%kappa
```

`predictors(i,j)` is observation `i` of macro factor `j`. An intercept is always
included internally, matching the upstream model requirement.

## Distribution example

```fortran
use vasicekfit, only : dp, vasicek_density, vasicek_cdf, vasicek_quantile

print *, vasicek_density(0.05_dp, 0.03_dp, 0.10_dp)
print *, vasicek_cdf(0.05_dp, 0.03_dp, 0.10_dp)
print *, vasicek_quantile(0.99_dp, 0.03_dp, 0.10_dp)
```

## Prediction conventions

`predict_link` returns the fitted probit-space linear predictor. `predict_response`
returns the conditional mean loss rate

```text
Phi(Phi^-1(p) + u' kappa).
```

`predict_quantiles` returns conditional Vasicek loss quantiles. By default the
results are on the loss-rate scale; pass `response_scale=.false.` for probit
space.

## HAC convention

The upstream R package delegates HAC long-run variance estimation to
`sandwich::lrvar`. This port supplies a dependency-free Bartlett/Newey-West
estimator. A lag can be supplied explicitly. Without one, the implementation
uses

```text
floor(4 * (N / 100)^(2 / 9)).
```

The statistical role is the same, but HAC results are not intended to be
bit-for-bit identical to every `sandwich` bandwidth, kernel, or prewhitening
configuration.

## Project layout

- `src/`: library modules
- `test/`: deterministic numerical tests
- `app/`: demonstration program
- `example/`: small usage examples
- `scripts/`: direct GNU Fortran validation scripts
- `original/`: unmodified upstream R package
- `provenance/`: supplied source archive and checksums

See `COVERAGE.md`, `PORTING_NOTES.md`, and `VALIDATION.md` for detailed mapping
and validation information.

## License

MIT. The original copyright and license are preserved. See `LICENSE` and
`NOTICE`.
