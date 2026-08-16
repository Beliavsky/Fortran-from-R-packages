# pracma-fortran

Modern Fortran/FPM translation of the computational core of the R package
`pracma` 2.4.6 (Practical Numerical Math Functions).

The library is organized as a numerical toolbox rather than a direct copy of
R's dynamic programming model.  It provides more than 400 public procedures
covering MATLAB-style array utilities, linear algebra, polynomials, special
functions, differentiation, integration, root finding, optimization,
interpolation, differential equations, signal/statistical routines, geometry,
and combinatorics.

## Build

```sh
fpm build
fpm test
fpm run --example example_linear_algebra
fpm run --example example_calculus
fpm run --example example_ode_signal
fpm run demo_pracma
```

GNU Fortran users can run the included validation scripts when FPM is not
available:

```sh
./scripts/test_gfortran.sh
./scripts/test_optimized.sh
```

On Windows, run `scripts\test_gfortran.bat` from a command prompt with
`gfortran` on `PATH`.

## Main modules

- `pracma_basic`: array creation/manipulation, MATLAB-style helpers, elementary
  functions, descriptive statistics, and coordinate-grid utilities.
- `pracma_linalg`: factorizations, linear solves, eigenanalysis, pseudoinverse,
  matrix functions, ranks, norms, and structured matrices.
- `pracma_polynomial`: polynomial arithmetic, roots, Chebyshev tools, Pade
  approximation, rational approximation, and companion matrices.
- `pracma_special`: gamma/beta-related functions, exponential/sine/cosine
  integrals, elliptic functions, zeta/polylogarithm, Lambert W, and related
  numerical functions.
- `pracma_differentiation`, `pracma_integration`, `pracma_roots`, and
  `pracma_optimization`: numerical calculus and solvers.
- `pracma_interpolation`: linear, nearest, spline, PCHIP, Akima, barycentric,
  and scattered/grid interpolation helpers.
- `pracma_ode`: fixed-step and adaptive ODE solvers plus structural and PDE
  time-stepping helpers.
- `pracma_signal_stats`: convolution, Fourier transforms, smoothing, peaks,
  outliers, entropy, Hurst estimates, and regression/statistical helpers.
- `pracma_geometry`: planar and spherical geometry, coordinate conversions,
  circle fitting, polygon/segment operations, and spatial sampling.
- `pracma_combinatorics`: permutations, combinations, number theory, random
  matrices, and matrix generators.
- `pracma_compat`: additional R/MATLAB-recognizable aliases and compatibility
  wrappers.

Use the aggregate module:

```fortran
use pracma
```

or import only the smaller module needed by an application.

## Scope and compatibility

The port represents 386 upstream function names directly after normalizing R
periods to Fortran underscores, and exposes 476 public Fortran names in total,
including typed callbacks, result types, status constants, overloads, and
compatibility helpers.  The compiled library deliberately omits plotting,
console/workspace manipulation, regular-expression/string utilities, and other
R-environment functions.  See `TRANSLATION_COVERAGE.md` and `PORTING.md` for
precise coverage and numerical adaptations.

The original R package sources, tests, manuals, data, demos, and metadata are
retained under `original/` for attribution and comparison.

## License

`pracma-fortran` is distributed under GPL-3.0-or-later, matching the upstream
`pracma` package.  The adapted quadratic-programming implementation originates
from GPL-2.0-or-later `quadprog` code and is therefore distributed here under
GPL-3.0-or-later as part of the combined work.  See `LICENSE` and `NOTICE.md`.
