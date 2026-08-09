# coneproj-fortran

Modern Fortran/FPM translation of the computational code in CRAN `coneproj` 1.23.

The original package implements primal/dual cone projections, a quadratic-programming transformation, constrained parametric regression, and shape-restricted regression. This translation exposes those numerical operations directly on Fortran arrays and omits R formula/model-frame/S3 printing and plotting infrastructure.

## Main API

```fortran
use coneproj

call cone_a(y, amat, ans)
call cone_b(y, delta, ans, vmat=vmat)
call qprog(q, c, amat, b, qp)
call make_delta(x, shape_increasing, delta, status)
call constreg_fit(y, xmat, amat, cr)
call shapereg_fit(y, x, shape_convex, sr)
call qr_decomp(xmat, qr)
call check_irreducible(edges, keep, reducible, equal_edges, status)
```

Constraint conventions match the R package:

- `cone_a`: rows of `amat` define the primal cone constraints.
- `cone_b`: columns of `delta` are cone edges; `vmat` spans the unconstrained linear subspace.
- `qprog`: minimizes `0.5*x'Q*x - c'x` subject to `amat*x >= b`.

Shape constants are `shape_increasing` through `shape_decreasing_concave` (1 through 8), matching the R helpers `incr`, `decr`, `conv`, `conc`, `incr.conv`, `decr.conv`, `incr.conc`, and `decr.conc`.

## Build

```text
fpm build
fpm test
```

No BLAS/LAPACK, R, Rcpp, or Armadillo dependency is required.

See `TRANSLATION_COVERAGE.md` and `VALIDATION.md` for implementation details and known differences.
