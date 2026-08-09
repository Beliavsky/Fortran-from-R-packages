# Dykstra-fortran

Modern Fortran translation of the computational core of the R package
`Dykstra` 1.0-0 by Nathaniel E. Helwig.

The package solves convex quadratic programs of the form

    minimize  0.5 * x^T D x - d^T x
    subject to A^T x >= b

with the first `meq` constraints optionally treated as equalities.  The solver
uses Richard L. Dykstra's cyclic projection algorithm, matching the upstream
R implementation rather than substituting another QP method.

## Features

- equality and inequality constraints in the upstream `solve.QP` orientation;
- positive-definite and positive-semidefinite quadratic matrices;
- diagonal fast path;
- source-compatible `factorized=.true.` input, where `Dmat` is `R^{-1}` and
  `R^T R = D`;
- source-compatible convergence tolerance rescaling and cycle count;
- self-contained symmetric Jacobi eigensolver; no BLAS/LAPACK dependency;
- explicit modern Fortran interfaces and an FPM package layout.

## Basic use

```fortran
use dykstra, only : dp, dykstra_result, dykstra_solve

type(dykstra_result) :: result

call dykstra_solve(dmat, dvec, amat, result, bvec=bvec, meq=1)
print *, result%solution
print *, result%value
```

`bvec`, `meq`, `factorized`, `maxit`, and `eps` are optional.  Omitted
`bvec` is treated as zero, just as in the R routine.

## Build

```text
fpm build
fpm test
```

The `scripts/` directory also contains strict GNU Fortran build scripts for
systems where FPM is unavailable.

## Translation scope

The numerical `dykstra()` routine is translated.  The R S3 print method and R
object/list infrastructure are not needed by the Fortran API.  See
`TRANSLATION_COVERAGE.md` for source-level fidelity details and documented
edge cases.

## License

The upstream package declares `GPL (>= 2)`.  The translated code is therefore
provided under GPL-2.0-or-later.  The upstream source is retained under
`original/Dykstra-master/`.
