# fitdistrplus-fortran

A modern Fortran 2018 and FPM computational port of the R package
`fitdistrplus` 1.2-6.

The library fits parametric distributions to uncensored or censored samples,
computes distribution diagnostics, and performs bootstrap refitting. Graphics,
R S3 objects, formulas, data frames, and expression capture are deliberately
outside the port.

## Implemented computations

- Maximum-likelihood estimation for ordinary samples
- Maximum-likelihood estimation for exact, left-, right-, and interval-censored observations
- Moment matching estimation (MME)
- Quantile matching estimation (QME), including weighted type-7 quantiles
- Maximum goodness-of-fit estimation (MGE): CvM, KS, AD, ADR, ADL, AD2R, AD2L, and AD2
- Maximum spacing estimation (MSE): KL, J, R, H, and V phi-divergences
- Automatic starting values for the built-in distributions
- Parameter-bound detection from distribution metadata
- Numerical observed-Hessian covariance for MLE fits
- Log likelihood, AIC, and BIC
- Descriptive mean, standard deviation, skewness, and kurtosis used by Cullen-Frey diagnostics
- KS, CvM, AD, and grouped chi-square goodness-of-fit statistics
- Parametric and nonparametric bootstrap refitting
- Censored-row bootstrap
- Bootstrap parameter quantiles and CDF confidence bands
- A self-contained Turnbull-style interval-censoring NPMLE

## Distribution interface

`type(distribution_model)` stores procedure pointers for log density, CDF,
quantile, raw moments, and random generation. User code can therefore add a
custom distribution without changing the fitting engine.

Built-in constructors are supplied for:

- normal and lognormal
- exponential, gamma, Weibull, uniform, logistic, Cauchy, and beta
- Poisson, geometric, and negative binomial (`size`, `mu`)

## Example

```fortran
program example
  use fitdistrplus
  implicit none
  type(distribution_model) :: dist
  type(fit_result) :: fit
  real(dp) :: x(6)

  x = [0.31_dp, 0.48_dp, 0.60_dp, 0.82_dp, 1.14_dp, 1.55_dp]
  call make_weibull(dist)
  call fitdist_auto(x, dist, method_mle, fit)
  write(*,'(a,2(1x,f10.5))') "shape, scale:", fit%estimate
end program example
```

## Build

With FPM:

```text
fpm test
fpm run demo_fitdistrplus
```

With GNU Make:

```text
make check
make optimized
make example
```

The checked target enables bounds, runtime, and backtrace checks. The optimized
target suppresses only two GNU Fortran descriptor-analysis warnings associated
with allocatable components of an `intent(out)` derived type.

## Numerical and interface differences

The R package delegates optimization and distribution definitions to the R
runtime. This port uses a self-contained bounded Nelder-Mead optimizer and a
Fortran callback model. It does not reproduce R named-list parameter fixing;
a fixed parameter can instead be embedded in a reduced custom callback.

Only MLE returns an observed-Hessian covariance matrix. The upstream asymptotic
MME covariance calculation is not reproduced. The interval-censoring NPMLE uses
the finite set of observed endpoints as support, which is deterministic and
self-contained but is not a line-for-line port of the survival-package engine.

See `PORTING_NOTES.md`, `API_MAP.md`, and `VALIDATION.md` for details.

## License

The upstream package declares `GPL (>= 2)`. This translation is distributed
under GPL-2.0-or-later. Upstream copyright notices and the original source
archive are retained under `upstream/`.
