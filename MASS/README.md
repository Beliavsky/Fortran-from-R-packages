# MASS modern Fortran

A self-contained modern Fortran 2018/FPM port of the computational portions of
MASS 7.3-66. The API replaces R formulas, S3 objects, and data frames with
explicit arrays and typed results.

## Implemented numerical areas

- Generalized inverses, null spaces, multivariate-normal simulation, 2-D KDE,
  numerical integration, rational approximation, contrasts, and contingency
  table conversion.
- Ordinary, generalized least-squares, ridge, robust M-estimation, LMS/LTS,
  residual diagnostics, Box-Cox/log transformations, and dose estimation.
- Classical and robust covariance, LDA, and QDA.
- Maximum-likelihood distribution fitting for normal, lognormal, Poisson,
  exponential, geometric, gamma, Weibull, beta, Cauchy, logistic,
  negative-binomial, Student-t, chi-square, and F distributions.
- Poisson and negative-binomial GLMs, log-linear models, negative-binomial
  simulation, and theta/gamma-shape estimation.
- Proportional-odds regression with logistic, probit, loglog, cloglog, and
  cauchit links.
- Correspondence analysis, multiple correspondence analysis, classical MDS,
  Sammon mapping, Kruskal nonmetric MDS, and Shepard disparities.
- UCV, BCV, and Sheather-Jones bandwidth selection.
- Numeric add/drop/stepwise AIC selection for supplied design-matrix columns.
- Negative-exponential self-start initialization.

## Deliberately omitted

Plotting and graphics helpers, datasets, formula parsing, S3 method dispatch,
R model frames, heterogeneous-list utilities, and `glmmPQL` are omitted.
`glmmPQL` depends on the full `nlme` mixed-effects object/runtime rather than a
standalone numerical kernel. `rms.curv` is also omitted because its public R
interface consumes nonlinear-model gradient and Hessian attributes managed by
R's `nls` object system.

The numerical `step_aic_linear`, `addterm_linear`, and `dropterm_linear`
routines operate on columns of an explicit matrix; they do not implement R's
formula hierarchy or interaction-scope rules.

## Build

With FPM:

```text
fpm test
fpm run
```

With GNU Make:

```text
make check
make optimized
make BUILD=build/check MODE=check example
```

GNU Fortran 14.2 was used for validation. The package does not require BLAS,
LAPACK, R, or C/C++ libraries.

## Minimal use

```fortran
use mass

type(regression_result) :: fit
call rlm_fit(x, y, fit, psi='bisquare')
```

See `app/demo_mass.f90`, `API_MAP.md`, and `PORTING_NOTES.md`.
