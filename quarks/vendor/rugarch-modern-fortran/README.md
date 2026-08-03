# rugarch-modern-fortran

**Status: experimental and incompletely validated.**

This is an independent modern Fortran 2018 translation of computational parts
of the R package **rugarch 1.5-6**. It is packaged as a standalone fpm library
and does not require R, Rcpp, Armadillo, or an external numerical library.
Plotting, R-specific S4/S3 infrastructure, parallel rolling, and resumable
rolling are intentionally omitted.

The code has passed strict GNU Fortran smoke tests and the included examples run,
but it has not been exhaustively checked against every `rugarch` function or edge
case. Verify results independently before production, scientific, financial, or
risk-management use.

## Implemented core

### Version 0.4 completed standalone numerical workflows

- Mean external regressors, a two-step ARCH-in-mean estimate, nonnegative
  variance-regressor effects, and variance-targeting controls through
  `fit_garch_extended` and `filter_garch_extended`.
- Automatic numerical Hessian, classical covariance, per-observation score,
  Newey-West sandwich covariance, and standard-error attachment for extended
  fits.
- Conditional raw-residual, kernel-smoothed, and semi-parametric-density
  bootstrap forecasts. Partial and full simulation/refit modes are available.
- Serial ARFIMA rolling forecasts, parametric simulation/refit distributions,
  residual bootstrap forecasts, multi-series fit/forecast workflows, and
  expanding-window cross-validation.
- Distribution skewness/kurtosis calculations, standardized GH parameter
  transformation, unconditional mean, volatility half-life, forecast
  performance statistics, and numeric forward-index utilities.

### Version 0.3 numerical workflows and diagnostics

- Numerical two-sided Hessians, classical covariance, Newey-West long-run
  covariance, and sandwich covariance from a supplied score matrix.
- Information criteria, weighted Ljung-Box/Monti tests, ARCH-LM, Nyblom
  parameter-stability statistics, sign-bias tests, and adjusted Pearson
  distribution goodness-of-fit tests.
- VaR duration, GMM moment, Hong-Li density, and Model Confidence Set tests.
- Squared, absolute, directional, trading-return, asymmetric, LINEX, and
  double-LINEX loss functions.
- Stationary and circular fixed-block bootstrap sampling.
- Multi-series fit/forecast wrappers, rolling one-step GARCH forecasts, and
  parametric simulation/refit distributions for fitted GARCH models.
- ARFIMA conditional forecasts and approximate AIC/BIC order selection.


- Standard GARCH, GJR-GARCH, eGARCH, APARCH, and integrated GARCH filtering,
  simulation, likelihood evaluation, fitting, and forecasting.
- Direct likelihood fitting for arbitrary requested orders of FIGARCH,
  component-GARCH, realGARCH, and Hentschel fGARCH. FIGARCH uses a configurable
  infinite-ARCH truncation.
- fGARCH submodels corresponding to rugarch's `GARCH`, `TGARCH`, `AVGARCH`,
  `NGARCH`, `NAGARCH`, `APARCH`, `ALLGARCH`, and `GJRGARCH` configurations.
- Joint realGARCH return and realized-measure likelihood, simulation, filtering,
  and measurement residuals.
- ARMA mean filtering inside volatility models.
- Persistence, unconditional variance, and news-impact calculations.
- Normal, skew-normal, standardized Student-t, skew Student-t, GED, skew GED,
  Johnson SU, standardized NIG, generalized hyperbolic, and generalized
  hyperbolic skew-t distributions.
- Density, CDF, quantile, random generation, and distribution maximum-likelihood
  fitting. GH-family CDFs and quantiles use native numerical quadrature.
- ARFIMA fractional differencing/integration, simulation, GPH estimation, and
  approximate AR fitting.
- VaR, expected shortfall, VaR coverage tests, ES test, directional-accuracy
  test, Berkowitz test, and quantile/VaR loss functions.

See [API_MAP.md](API_MAP.md) for exact coverage and omissions.

## Build and run

```console
fpm build
fpm test
fpm run
fpm run --example fit_csv -- data/sample_returns.csv
fpm run --example complete_workflows
```

On Windows, `run_examples.bat` runs the demonstration, CSV example, and completed-workflows example.

Plain `fpm run` runs only the no-argument demonstration. The CSV program is an
fpm example rather than a normal executable, so it is not launched accidentally.

## Basic GARCH use

```fortran
program example
   use rugarch
   implicit none

   type(garch_fit_result) :: fit
   real(dp) :: returns(1000)

   ! Fill returns first.
   fit = fit_garch11(returns,dist_std,fit_shape=.true.)

   print *,fit%spec%omega
   print *,fit%spec%alpha
   print *,fit%spec%beta
end program example
```

## Advanced model fitting

```fortran
! FIGARCH(2,d,1), using 500 infinite-ARCH coefficients.
fit = fit_figarch(returns,p=2,q=1,truncation=500)

! Component GARCH(1,1).
fit = fit_csgarch11(returns)

! realGARCH(1,1); realized must be positive and match returns in length.
fit = fit_realgarch11(returns,realized)

! Hentschel ALLGARCH(1,1). Other fgarch_* constants select restrictions.
fit = fit_fgarch11(returns,submodel=fgarch_allgarch)
```

## Completed workflow API

```fortran
! External mean/variance regressors, variance targeting, and covariance output.
extended = fit_garch_extended(returns,model_sgarch,1,1,dist_std, &
   mean_regressors=mx,variance_regressors=vx,variance_targeting=.true., &
   fit_shape=.true.)

! Conditional bootstrap forecast density.
boot = garch_bootstrap_forecast(extended%fit,10,500,sampling_kernel, &
   bootstrap_partial)

! Serial ARFIMA rolling and parameter-distribution workflows.
roll = rolling_arfima_forecast(returns,1,1,250,refit_every=25)
dist = arfima_parametric_distribution(arfima_fit,500)
```

## Model constants

- `model_sgarch`, `model_gjrgarch`, `model_egarch`
- `model_aparch`, `model_igarch`, `model_figarch`
- `model_csgarch`, `model_realgarch`, `model_fgarch`

fGARCH submodel constants are:

- `fgarch_garch`, `fgarch_tgarch`, `fgarch_avgarch`
- `fgarch_ngarch`, `fgarch_nagarch`, `fgarch_aparch`
- `fgarch_allgarch`, `fgarch_gjrgarch`

## Distribution constants

- `dist_norm`, `dist_snorm`
- `dist_std`, `dist_sstd`
- `dist_ged`, `dist_sged`
- `dist_jsu`
- `dist_nig`, `dist_ghyp`, `dist_ghst`

For `dist_ghyp`, `garch_spec%lambda` is the generalized-hyperbolic lambda
parameter. This is separate from `garch_spec%fgarch_lambda`, the Hentschel
fGARCH variance power.

## Meaning of "full fitting"

For FIGARCH, component-GARCH, realGARCH, and fGARCH, the project now estimates
the model-specific parameters used by the implemented recursions rather than
providing filter-only placeholders. It supports arbitrary positive `p` and `q`
orders and all listed fGARCH submodel restrictions.

The dependency-free API does not reproduce R's object-oriented fitting
environment or every solver and initialization control. The standalone
`fit_garch_extended` workflow adds external regressors, ARCH-in-mean,
variance targeting, and automatic classical/robust covariance output. Its
regressor treatment is an experimental numerical equivalent rather than a
bit-for-bit reproduction of R's joint solver stack.

## Limitations

- The dependency-free Nelder-Mead optimizer is less sophisticated than the
  solver stack used by R `rugarch`; convergence status and estimates must be
  inspected.
- External-regressor and ARCH-in-mean estimation use a stable two-step
  regression/volatility procedure. Variance regressors are estimated
  nonnegatively and then applied to conditional variance. This is not identical
  to every R `rugarch` joint-likelihood configuration.
- Numerical covariance matrices can be singular or poorly conditioned near
  parameter boundaries. `covariance_status` must be checked.
- Raw, kernel, and semi-parametric bootstrap forecasts are implemented. The
  semi-parametric tail sampler is dependency-free and therefore differs from
  the R `spd` package implementation.
- GH-family CDFs, quantiles, and numerical moments are slower than normal, t,
  or GED functions.
- FIGARCH is a finite-truncation implementation of the infinite ARCH form.
- Serial rolling and multi-series workflows are implemented. Parallel rolling
  and checkpoint/resume support are intentionally omitted.
- R S4/S3 classes, methods, reports, and plotting are intentionally omitted.
- Numerical equivalence to all R `rugarch` combinations is not established.

## License and attribution

The original R package is GPL-3. This translation is distributed under
GPL-3.0-only. Original authorship and provenance are retained in `NOTICE`,
`ORIGIN.md`, and `reference/`.
