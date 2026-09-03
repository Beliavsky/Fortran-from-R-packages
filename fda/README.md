# fda — modern Fortran computational translation

This directory is a modern free-form Fortran/FPM translation of a portable
computational subset of the R package **fda 6.3.0** (2025-05-21), authored by
James Ramsay with contributions from Giles Hooker and Spencer Graves.

The upstream package supports the functional-data methods described by Ramsay
and Silverman (2005) and Ramsay, Hooker, and Graves (2009).  This translation
focuses on numerical algorithms that have a natural Fortran API and omits R
S3/formula/plotting/data-set interface machinery.

## Implemented numerical areas

- constant, B-spline, Fourier, monomial, exponential, power, and polygonal
  bases, including integer derivatives;
- B-spline endpoint multiplicities compatible with `bsplineS` and numerical
  derivative/Gram/roughness matrices;
- functional-data coefficient objects, evaluation, means, centering, and
  cross-function inner products;
- penalized basis smoothing with vector observation weights, fitted values,
  SSE, effective degrees of freedom, GCV, and the observation-to-coefficient
  map;
- `lambda2df`, `df2lambda`, `lambda2gcv`, and `project.basis` equivalents;
- regularized univariate functional PCA and functional CCA;
- generalized metric SVD (`geigen`), symmetric solves, trapezoidal matrix
  integration, Simpson quadrature, polynomial extrapolation, and zero-bracketing;
- the package's homogeneous linear differential-equation conversion and
  adaptive Cash–Karp solver (`derivs`, `odesolv`, `rkqs`, `rkck`).

See `API_COVERAGE.md` for exact mappings and the remaining parity targets.

## Dependencies

Place this directory beside the shared packages in
`Fortran-from-R-packages`:

```text
Fortran-from-R-packages/
├── fda/
├── rfortran-core/
└── rfortran-linalg/
```

The FPM manifest uses only sibling path dependencies:

```toml
rfortran-core   = { path = "../rfortran-core" }
rfortran-linalg = { path = "../rfortran-linalg" }
```

`rfortran-core` supplies the shared `dp` kind.  `rfortran-linalg` supplies
linear systems, Cholesky factors, symmetric eigenanalysis, inverses, and SVD;
that shared package in turn pins `fortran-lapack`, so this translation neither
vendors BLAS/LAPACK nor links system BLAS/LAPACK.

Although upstream `DESCRIPTION` imports `fds` and `deSolve`, the translated
portable kernels do not require either namespace.  No `fds` call occurs in the
upstream R computational files used here, and `odesolv.R` contains its own
adaptive Runge–Kutta implementation, which is translated directly.

## Build and test

From this directory, with the sibling dependencies present:

```text
fpm build
fpm test
fpm run --example smooth_sine
```

The intended compiler is GNU Fortran on Windows or another standards-conforming
Fortran 2018 compiler.  No separately installed BLAS or LAPACK library should
be required because `rfortran-linalg` uses its pinned pure-FPM dependency.

## Example

`example/smooth_sine.f90` builds a cubic B-spline basis on `[0,1]`, smooths a
sampled sine curve with a second-derivative penalty, and prints the effective
degrees of freedom and SSE.

## Numerical note on `rkck`

Upstream `R/odesolv.R` has `575/512` multiplying the second Cash–Karp stage in
the sixth-stage state.  The Cash–Karp tableau requires `175/512`.  The Fortran
translation intentionally uses `175/512`: retaining the upstream literal makes
the adaptive solver shrink steps excessively and materially degrades a simple
`u'' + u = 0` reference solution.  This deliberate correction is also recorded
in `NOTICE.md` and `API_COVERAGE.md`.

## License and provenance

Upstream declares `GPL (>= 2)`.  This translation is therefore distributed
under **GPL-2.0-or-later**.  `COPYING` and `LICENSE` contain GPL version 2;
`upstream/` preserves the package `DESCRIPTION`, `NAMESPACE`, `MD5`, and NEWS
metadata from the supplied source archive.  See `NOTICE.md` for detailed
provenance and dependency decisions.
