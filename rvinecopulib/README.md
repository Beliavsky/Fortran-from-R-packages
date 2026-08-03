# rvinecopulib-fortran

A self-contained modern Fortran/FPM translation of the continuous parametric
modeling core of `rvinecopulib`.

The library provides bivariate copula models and fitted C-vine and D-vine
models without requiring R, Rcpp, Eigen, Boost, or the C++ `vinecopulib`
backend.

## Implemented pair-copula families

- independence
- Gaussian
- Student-t
- Clayton
- Gumbel
- Frank
- Joe
- BB1, BB6, BB7, and BB8
- Tawn
- rotations 0, 90, 180, and 270 degrees

Each `bicop_model` supports density, distribution, the two h-functions,
numerical inverse h-functions, simulation, Kendall's tau, log likelihood,
AIC, and BIC. `fit_bicop` performs bounded maximum likelihood estimation and
`select_bicop` performs AIC/BIC family selection.

## Vine models

`cvine_model` and `dvine_model` support:

- arbitrary variable order
- triangular pair-copula stores
- sequential maximum-likelihood fitting
- pair-family selection
- density and log density
- Rosenblatt and inverse Rosenblatt transforms
- simulation
- Monte Carlo CDF evaluation
- truncation
- AIC and BIC

Fortran data matrices use **variables by observations**. Thus a bivariate
sample has shape `(2, n)` and a `d`-dimensional sample has shape `(d, n)`.

## Build

```text
fpm build
fpm test
fpm run
```

The version in `fpm.toml` is the FPM-compatible numeric value `0.7.3`; the
attached upstream package identifies itself as 0.7.3.1.0.

GNU Fortran scripts are also supplied:

```text
scripts/build_all.sh checked
scripts/build_all.sh optimized
```

On Windows with GNU Fortran:

```text
scripts\build_all.bat checked
scripts\build_all.bat optimized
```

## Minimal example

```fortran
program example
  use rvinecopulib
  implicit none
  type(bicop_model) :: copula

  copula = make_bicop(bicop_clayton, 0, [2.0_dp])
  print *, copula%pdf(0.3_dp, 0.7_dp)
  print *, copula%cdf(0.3_dp, 0.7_dp)
  print *, copula%tau()
end program
```

See `example/` and `app/demo_rvinecopulib.f90` for fitting and vine examples.

## Scope differences

This port intentionally omits R class methods, plotting, formula/data-frame
handling, nonparametric TLL copulas, mixed discrete/continuous likelihoods,
arbitrary regular-vine matrix validation and Dissmann structure selection,
quasi-random sequences, parallel execution, JSON serialization, and fitted
marginal-distribution objects. C-vines and D-vines are implemented natively and
cover the most commonly used structured-vine workflows.

See `API_MAP.md` and `PORTING_NOTES.md` for details.
