# limSolve-fortran

Modern Fortran/FPM translation of the computational core of the R package
`limSolve` 2.0.3.

The public API works directly with numeric arrays.  R model objects, names,
warnings, printing, and plotting are intentionally not reproduced.

## Implemented computational API

- generalized-inverse solution of `A X = B` (`Solve` -> `solve_generalized`)
- nonnegative least squares (`nnls`)
- least-distance programming (`ldp`)
- equality/inequality least-distance problems (`ldei`)
- constrained least squares (`lsei`)
- linear programming (`linp`)
- equation/variable resolution diagnostics (`resolution`)
- tridiagonal, banded, and almost-block-diagonal linear systems
- `xranges`, `varranges`, and `varsample`
- `xsample` with mirror, random-direction, and coordinate-direction walks
- lower and upper bounds expressed in the same `G*x >= H` convention

## Dependencies

The original R package imports `lpSolve` and `quadprog`.  This release vendors
and links the Fortran translations produced in the same translation project:

```toml
[dependencies]
lpsolve-fortran = { path = "vendor/lpSolve-fortran" }
quadprog-fortran = { path = "vendor/quadprog-fortran" }
```

`linp` and range calculations use the translated `lpSolve` simplex/MILP core.
The convex QP step used by `ldp` and `lsei` uses the translated Goldfarb-Idnani
`quadprog` core.

## Build

```text
fpm build
fpm test
```

A strict GNU Fortran build is available in `scripts/test_gfortran.sh` and
`scripts/test_gfortran.bat`.

See `TRANSLATION_COVERAGE.md` for the numerical/architectural differences from
the R package and `UPSTREAM_PROVENANCE.md` for source provenance.
