# Translation coverage

## Summary

The upstream package contains 498 unique top-level R function definitions
identified from its source files.  Many are internal callbacks, plotting
helpers, console/workspace commands, regular-expression utilities, or dynamic
closure factories rather than reusable numerical algorithms.

The Fortran library exposes 476 public names.  After case-folding and mapping R
periods/Fortran underscores, 386 upstream names have direct public name
matches.  Additional public names are typed callback interfaces, result types,
status constants, overload-specific procedures, and Fortran compatibility
helpers.

## Numerical families translated

| Family | Coverage |
|---|---|
| MATLAB-style vectors, matrices, grids, reshaping, indexing | Broad |
| Descriptive statistics and elementary numerical helpers | Broad |
| Dense linear algebra, matrix decompositions/functions | Broad |
| Polynomial and rational approximation | Broad |
| Special mathematical functions | Broad |
| Numerical differentiation and Jacobians/Hessians | Broad |
| One-, two-, and three-dimensional quadrature | Broad |
| Scalar and multivariate root finding | Broad |
| Unconstrained, bounded, constrained, LSQ, LP, NNLS, and QP optimization | Broad |
| Linear/spline/PCHIP/Akima/barycentric interpolation | Broad |
| ODE, structural dynamics, and parabolic PDE stepping | Broad, with shared RK backend for some aliases |
| Signal processing and time-series/statistical helpers | Broad |
| Planar/spherical geometry and spatial sampling | Broad |
| Number theory, combinatorics, sorting, and matrix generators | Broad |

## Not compiled

The following categories remain in `original/` but are not represented as
compiled numerical procedures:

- graphics and plotting;
- R workspace/environment/filesystem/console operations;
- regular-expression and string-manipulation helpers;
- formula/data-frame/S3 behavior;
- functions whose principal result is an R closure;
- small upstream-local helper names that have no standalone API meaning.

## Adapted rather than line-identical

- FFT uses a portable direct DFT.
- Several named adaptive ODE methods share RKF45.
- SVD/pseudoinverse and some matrix functions use symmetric eigendecomposition.
- Optimization methods use self-contained portable algorithms.
- Random streams are reproducible within this library but do not match R.
- The generalized R shooting/BVP framework is represented by a scalar
  finite-difference BVP solver.

See `PORTING.md` for details and `API.md` for the public procedure inventory.
