# NOTICE

This directory is a modern Fortran translation of computational code from the
R package `randompack` 0.1.10 by Kristjan Jonasson, University of Iceland.
The pristine upstream source used for the translation is retained under
`upstream/`, including its `DESCRIPTION`, `LICENSE`, source comments, and
third-party notices. The upstream `inst/THIRD-PARTY-NOTICES` file is also
copied verbatim to the package root as `THIRD-PARTY-NOTICES`.

The translation retains the upstream BSD-3-Clause licensing terms. Some RNG
algorithms originate from separately credited public-domain, BSD, 0BSD/MIT-0,
MIT, CC0, and other sources listed in `THIRD-PARTY-NOTICES`. Those notices and
attributions remain applicable to translated or adapted algorithmic material.

The Fortran package uses the sibling `rfortran-linalg` package for the ordinary
Cholesky operation needed by multivariate-normal simulation. The maintained
Fortran source also translates the upstream package's own `rp_dpstrf.c`
complete-pivoting positive-semidefinite Cholesky fallback. It does not copy
BLAS, LAPACK, `rfortran-linalg`, or another translated R package into this
directory. The linear-algebra backend is owned by the sibling
`rfortran-linalg` package; `randompack` neither selects nor vendors a
BLAS/LAPACK implementation directly.

Architecture-specific AVX2, AVX-512, NEON, SLEEF, and R interface code is not
copied into the maintained Fortran implementation. The three upstream SIMD RNG
engines are reproduced portably by scalar emulation of their eight logical
streams, preserving stream semantics rather than machine-vector acceleration.

The maintained `src/randompack_openlibm.f90` translates the double-precision routines retained in upstream `src/openlibm.inc`, which identifies its OpenLibm/fdlibm-derived code as BSD-2-Clause. See `THIRD-PARTY-NOTICES`.
