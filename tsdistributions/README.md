# tsdistributions-fortran

Modern Fortran 2018 computational translation of R package `tsdistributions` 1.0.4.

## Scope

The library implements standardized Normal, Student-t, Fernandez-Steam skew-normal and skew-t, GED and skew-GED, NIG, generalized hyperbolic, Johnson SU, and generalized-hyperbolic skew-t distributions. It includes densities, CDFs, quantiles, random generation, theoretical skewness and kurtosis, constrained moment-domain calculations, numerical maximum-likelihood fitting, Hessian/OPG/sandwich covariance estimates, simulation profiles, and the package's semiparametric GPD-tail/kernel-interior model.

R/TMB classes, data tables, plotting, formula dispatch, parallel futures, and automatic-differentiation runtime integration are intentionally omitted. Typed derived types and explicit procedures replace R objects.

## Build

```sh
make checked
make optimized
```

With FPM installed:

```sh
fpm test
fpm run --example example_tsdistributions
```

## Minimal use

```fortran
use tsdistributions

type(distribution_parameters) :: p
real(dp) :: density

p = distribution_parameters(mu=0.0_dp, sigma=1.0_dp, shape=8.0_dp)
density = ddist('std', 0.0_dp, p)
```

See `example/example_tsdistributions.f90`, `API_MAP.md`, and `PORTING_NOTES.md`.
