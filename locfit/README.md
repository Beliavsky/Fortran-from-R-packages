# locfit-fortran

`locfit-fortran` is a modern free-format Fortran/FPM translation of the
computational core of the R package **locfit 1.5-9.12** by Catherine Loader
and contributors.

The upstream package implements local regression, local likelihood and local
density estimation methods described in Loader's work.  This translation
focuses on numerical/statistical functionality and deliberately omits R
formula/S3 infrastructure and plotting/lattice code.

## License and provenance

The upstream `locfit` package declares `GPL (>= 2)` and contains Lucent
Technologies copyright notices.  The upstream README records the package's
2005 relicensing.  This translation is therefore distributed under
GPL-2.0-or-later.  See `COPYING` and `LICENSE_NOTICE`.

A complete unmodified copy of the upstream package used for the translation
is retained in `upstream/locfit-R/`.

## Implemented in v0.1.0

- Local polynomial regression and local likelihood at arbitrary evaluation
  points.
- Spherical and product kernels, nearest-neighbor plus fixed bandwidths,
  predictor styles, and automatic scale selection.
- Polynomial fitting bases through the degrees supported by upstream locfit's
  standard basis code (full spherical basis through degree 3; product basis
  through the configured degree).
- Rectangular, Epanechnikov, biweight, tricube, triweight, Gaussian,
  triangular, quartic, six-cube, exponential, Maclean and parametric kernels.
- Gaussian, binomial/logistic, robust-binomial, Poisson, Gamma, geometric,
  Weibull, circular, Huber-robust and Cauchy local-likelihood family terms,
  with default/canonical links and the censored likelihood branches present in
  upstream `family.c`.
- Local coefficient covariance matrices and pointwise standard errors.
- Derivative evaluation from the fitted local polynomial.
- One-dimensional kernel density estimation and one-dimensional local
  log-density estimation.
- KDE bandwidth criteria/selectors: AIC, ordinary CV, LSCV, BCV,
  Sheather-Jones plug-in and GKK.
- Linear, cubic Hermite and rectangular-cell interpolation primitives.
- Raw/deviance/Pearson/likelihood residual helpers, studentization helper,
  GCV/AIC/Cp formulas, Kaplan-Meier mean residual life.
- High-level iterative robust and quasi-likelihood fitting helpers corresponding
  to the computational loops in `locfit.robust` and `locfit.quasi`.

See `docs/SOURCE_MAP.md` and `docs/TRANSLATION_STATUS.md` for exact coverage.

## Build with FPM

```text
fpm build
fpm test
fpm run --example smooth_example
```

The package has no external Fortran dependencies.

## Minimal example

```fortran
program demo
  use locfit
  implicit none

  integer, parameter :: n=21, m=9
  real(dp) :: x(n,1), y(n), xe(m,1)
  type(locfit_options) :: opt
  type(locfit_result) :: fit
  integer :: i

  do i=1,n
    x(i,1)=-2.0_dp+4.0_dp*real(i-1,dp)/real(n-1,dp)
    y(i)=sin(x(i,1))
  end do
  do i=1,m
    xe(i,1)=-2.0_dp+4.0_dp*real(i-1,dp)/real(m-1,dp)
  end do

  opt%nn=0.55_dp
  opt%degree=2
  opt%kernel=wtcub
  call locfit_fit(x,y,xe,fit,opt)

  do i=1,m
    write(*,'(3f14.7,i4)') xe(i,1),fit%fit(i),fit%se(i),fit%status(i)
  end do
end program demo
```

## Main public API

The convenience module `locfit` re-exports the package modules.  The primary
entry points are:

- `locfit_fit`
- `locfit_fit_one`
- `locfit_derivative_at`
- `locfit_robust_fit`
- `locfit_quasi_fit`
- `kernel_density_1d`
- `local_density_1d`
- `kde_select`
- `kde_criterion`

Configuration is supplied using `type(locfit_options)`; results are returned
in `type(locfit_result)`.

## Validation

`test/test_locfit.f90` exercises:

- kernel and special-function identities,
- exact recovery of univariate and multivariate polynomials,
- derivative recovery,
- binomial and Poisson intercept-only likelihood fits,
- equality of degree-zero local log-density and the corresponding KDE,
- bandwidth criteria/selectors,
- Hermite interpolation,
- residual/criterion helpers,
- automatic predictor scaling, and
- robust outlier downweighting.

The translation was also compiled and run with gfortran using bounds/runtime
checking, floating-point traps, warnings, and implicit-interface errors.  FPM
was not installed in the translation environment, so the FPM layout was
validated by compiling the same source/test units directly with gfortran.
