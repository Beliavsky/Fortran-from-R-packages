# fcopulae-modern-fortran

A modern Fortran computational translation of the R package `fCopulae`.

The project translates the package's numerical copula routines into plain-array Fortran APIs. Plotting, sliders, S4 objects, `timeSeries` metadata, and R-specific formatting are intentionally excluded.

## Implemented numerical scope

### Archimedean copulas

- All 22 generator families listed by `archmList()`
- Default parameters, parameter ranges, and validation
- Generator, inverse generator, endpoint value, and numerical first/second derivatives
- Kendall distribution function and inverse
- CDF, density, and random generation
- Kendall's tau, Spearman's rho, and lower/upper tail coefficients
- Bounded maximum-likelihood fitting with AIC and BIC

### Elliptical copulas

- Normal, Cauchy, Student-t, logistic, Laplace, Kotz, and exponential-power families
- Elliptical generators and normalization constants
- Univariate marginal PDF, CDF, and quantile functions
- Bivariate joint densities and copula CDFs/densities
- Random generation for every family
- Kendall's tau, Spearman's rho, and tail coefficients
- Bounded maximum-likelihood fitting, including Student degrees of freedom and Kotz/exponential-power shape parameters

### Extreme-value copulas

- Gumbel
- Galambos
- Husler-Reiss
- Tawn
- BB5
- Pickands dependence functions and numerical derivatives
- CDF, density, conditional-inversion simulation
- Kendall's tau, Spearman's rho, and tail coefficients
- Bounded maximum-likelihood fitting
- Pickands/CDF evaluation for the package's supplemental Gumbel-II, independence, and comonotonic cases

### Empirical and limiting copulas

- Empirical copula CDF at arbitrary points
- Empirical copula grids and finite-difference density grids
- Frechet upper, independence, lower, and PSP copulas
- Marshall-Olkin CDF
- Debye functions of positive integer order
- Pseudo-observations and sample Kendall/Spearman statistics

## Build

GNU Fortran, LAPACK, and BLAS are required.

```sh
make debug
make release
```

The debug target uses Fortran 2018, warnings as errors, bounds/runtime checking, and backtraces. The release target uses `-O2` with the same warning policy.

An `fpm.toml` file is included. It was not used for the recorded validation because `fpm` was unavailable in the validation environment.

## Applications

```sh
build/bin/demo_fcopulae
build/bin/dependence_example
build/bin/fit_csv data/sample_uv.csv archm 1 500
build/bin/fit_csv data/sample_uv.csv elliptical t 500
build/bin/fit_csv data/sample_uv.csv ev gumbel 500
```

The CSV file may contain either `U,V` or `Date,U,V`. Both pseudo-observations must lie strictly inside `(0,1)`.

## Numerical implementation notes

- Bounded Nelder-Mead replaces R's `nlminb` orchestration.
- Numerical quadrature and finite differences replace several R `integrate`, interpolation, and derivative paths.
- Normal, Cauchy, and Student-t elliptical CDFs use dedicated bivariate routines inherited from the translated `fMultivar` numerical layer.
- Logistic, Laplace, Kotz, and exponential-power CDFs use adaptive numerical integration.
- For the four nonstandard elliptical families, simulation draws the correct radial elliptical latent distribution and converts each margin to pseudo-uniform ranks. This preserves the copula sample ordering but is not an exact analytic marginal transform.
- The generic Archimedean random generator uses the Kendall-distribution construction for every family rather than retaining separate R-only fast paths.

## Source corrections

Two apparent inconsistencies in the attached R source were not copied blindly:

1. The documented Tawn parameter range requires both asymmetry parameters in `[0,1]`, but the original default sets the first one to `2`. The Fortran default is `[0.5, 0.5, 2]`.
2. The original internal Marshall-Olkin routine uses `max` between two branches, which can exceed `min(u,v)` and therefore violate a copula bound. The Fortran implementation uses the standard `min` formula. The original `.dmoCopula` duplicates that probability expression rather than computing a density, so no misleading continuous-density routine is exposed for this copula with a singular component.

## Licensing

The original package declares `GPL (>= 2)`. This translation uses `SPDX-License-Identifier: GPL-2.0-or-later` in every Fortran source, application, example, and test file. `LICENSE` contains GNU GPL version 2.

## Explicit exclusions

- Plotting, contour/perspective helpers, sliders, and GUI code
- S4 `fCOPULA` objects and R show/print infrastructure
- R list/grid adapters and attributes
- `timeSeries`/date metadata
- Exact R optimizer, quadrature, interpolation, and random-stream equivalence
- Separate exposure of alternative direct-formula wrappers when the generic generator/Pickands formula computes the same copula
- Parallel fitting or vectorized R callback infrastructure

See `API_MAP.md` and `VALIDATION.md` for routine-level mapping and executed tests.
