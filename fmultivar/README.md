# fmultivar-modern-fortran

A modern Fortran translation of the computational routines in the R package
`fMultivar` version 4031.84.

The translation excludes plotting and R class infrastructure. Numerical results
are returned as scalars, allocatable arrays, or plain Fortran derived types.

## License

The original package declares `GPL (>= 2)`. This project therefore uses
`SPDX-License-Identifier: GPL-2.0-or-later`. The complete GNU GPL version 2 text
is in `LICENSE`. Every Fortran source file contains the SPDX identifier and an
explicit version 2-or-later notice. `make check` verifies the headers.

## Implemented computational features

### Bivariate distributions

- Bivariate Normal density, CDF, and simulation
- Bivariate Student-t density, CDF, and simulation
- Bivariate Cauchy density, CDF, and simulation
- Elliptical densities for:
  - Normal
  - Cauchy
  - Student-t
  - Logistic
  - Laplace
  - Kotz
  - Exponential power

The Normal and Student-t CDFs use adaptive conditional one-dimensional
quadrature. The exact origin identity for the bivariate Normal CDF is handled
analytically.

### Multivariate distributions

- Multivariate Normal density, log density, simulation, rectangular
  probability, and equicoordinate quantile
- Multivariate Student-t density, log density, simulation, rectangular
  probability, and equicoordinate quantile
- Multivariate skew-Normal density, simulation, and rectangular probability
- Multivariate skew-Student density, simulation, and rectangular probability
- Multivariate skew-Cauchy density, simulation, and rectangular probability
- Compatibility procedures corresponding to the original names:
  `dmvnorm`, `pmvnorm`, `qmvnorm`, `rmvnorm`, `dmvt`, `pmvt`, `qmvt`, `rmvt`,
  `dmsn`, `pmsn`, `rmsn`, `dmst`, `pmst`, `rmst`, `dmsc`, `pmsc`, and `rmsc`
- Historical wrapper aliases `dmvsnorm`, `pmvsnorm`, `rmvsnorm`, `dmvst`,
  `pmvst`, and `rmvst`

Higher-dimensional rectangular probabilities are estimated by reproducible
Monte Carlo simulation. Probability procedures return both the estimate and
its Monte Carlo standard error. Equicoordinate quantiles use common random
numbers during bisection.

### Distribution fitting

- Multivariate Normal sample mean and covariance fitting
- Skew-Normal maximum-likelihood fitting
- Skew-Student maximum-likelihood fitting with fixed or estimated degrees of
  freedom
- Skew-Cauchy maximum-likelihood fitting
- Direct parameterization with a positive-definite scale matrix constructed
  from an unconstrained Cholesky factor
- Bounded Nelder-Mead optimization
- Numerical Hessian and parameter covariance matrix
- Compatibility wrappers `msn_fit`, `mst_fit`, `msc_fit`, and `mvfit`

The fit result contains location, scale matrix, shape vector, degrees of
freedom, log likelihood, numerical Hessian, covariance matrix, iteration count,
and convergence flag.

### Integration and data utilities

- Adaptive Simpson one-dimensional integration
- The original nine-point `integrate2d` unit-square rule, exposed as both
  `integrate2d_rule` and `integrate2d`
- Adaptive Gauss-Legendre two-dimensional integration over arbitrary finite
  rectangles
- Halton quasi-Monte Carlo integration in 1 through 20 dimensions, exposed as
  `adapt_integrate_nd` and the compatibility wrapper `adapt`
- Cartesian grid construction
- Plain `grid_data` construction through `make_grid_data` and `griddata`
- Axis-aligned bivariate Normal kernel density estimation
- Two-dimensional histograms
- Square binning with cell centers and centers of mass through
  `square_binning` and `squarebinning`
- Hexagonal binning with cell centers and centers of mass through
  `hex_binning` and `hexbinning`

The square-binning implementation uses explicit `bins_x` by `bins_y` cells and
retains observations on the maximum boundaries. This fixes the edge-loss risk
in the historical R indexing while preserving the intended calculation.

### Numerical support

- Reproducible uniform, Normal, Gamma, chi-squared, multivariate Normal, and
  multivariate Student-t random generation
- Normal and Student-t PDF, CDF, and quantile functions
- Regularized incomplete beta evaluation
- LAPACK-backed Cholesky factorization, positive-definite solves and inverses,
  general matrix inversion, log determinants, and sample covariance
- A reusable bounded Nelder-Mead optimizer

## Build

GNU Fortran, LAPACK, and BLAS are required.

```sh
make check
make release-check
```

The debug workflow uses:

```text
-std=f2018 -O0 -g -Wall -Wextra -Wimplicit-interface -Werror
-fcheck=all -fbacktrace
```

The optimized workflow uses:

```text
-std=f2018 -O2 -Wall -Wextra -Wimplicit-interface -Werror -fbacktrace
```

An `fpm.toml` manifest is included and links LAPACK and BLAS. `fpm` was not
available in the validation environment, so the manifest is not claimed as
tested.

## Applications

Run the demonstration after `make check`:

```sh
build/debug/demo_fmultivar
```

Fit a dated CSV file whose first column is a date and remaining columns are
numeric observations:

```sh
build/debug/fit_csv data/sample_returns.csv normal
build/debug/fit_csv data/sample_returns.csv snorm 0 900
build/debug/fit_csv data/sample_returns.csv st 6 900
build/debug/fit_csv data/sample_returns.csv st 0 1200
build/debug/fit_csv data/sample_returns.csv cauchy 0 900
```

For `st`, a positive third argument fixes the degrees of freedom. A nonpositive
value estimates it.

## Important numerical differences from R

- R's `mvtnorm` and `sn` probability engines are not reproduced exactly.
  Higher-dimensional and skewed rectangular probabilities use reproducible
  Monte Carlo simulation and report simulation error.
- R's `sn::msn.mle` and `sn::mst.mple` are replaced by bounded Nelder-Mead
  likelihood optimization with finite-difference Hessians. Exact optimizer
  endpoints, penalties, and standard errors are not claimed to match R.
- `cubature::adaptIntegrate` is represented by adaptive 2-D Gauss-Legendre
  integration and a Halton quasi-Monte Carlo routine for general dimensions.
  Vector-valued integrands and the full `cubature` option set are not included.
- Random streams do not reproduce R's RNG.
- Functions require finite, already-cleaned numeric arrays. Automatic R-style
  coercion, recycling, missing-value handling, and attributes are not included.

## Explicit exclusions

- Plot, contour, perspective, image, rug, and slider functions
- S3 and S4 classes and methods, including `fDISTFIT`, `gridData`,
  `squareBinning`, and `hexBinning` class machinery
- R formulas, calls, titles, descriptions, and package metadata objects
- Exact `mvtnorm`, `sn`, and `cubature` algorithms
- R graphics and interactive Tcl/Tk infrastructure

See `API_MAP.md` and `VALIDATION.md` for exact coverage and test details.
