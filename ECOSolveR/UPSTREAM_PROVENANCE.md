# Upstream provenance

## Supplied package

The user supplied `ECOSolveR-master.zip`, containing ECOSolveR version 0.6.1 source, R wrappers, C interface code, tests, and test data.

The package declares:

```text
License: GPL(>=3)
```

The original supplied tree is preserved under `original/ECOSolveR-master/`.

## ECOS solver submodule

The supplied archive's `.gitmodules` contains:

```text
[submodule "src/ecos"]
    path = src/ecos
    url = https://github.com/bnaras/ecos.git
    branch = r-patches
```

but the uploaded ZIP contains an empty `src/ecos/` directory.  The referenced `r-patches` branch was therefore consulted for ECOS's cone definitions, solver/KKT architecture, ECOS_BB behavior, and licensing.

## Sparse LDL provenance

Upstream ECOS uses the sparse LDL package by Timothy A. Davis.  The v0.2-v0.4 Fortran sparse backend follows the same published LDL symbolic/numeric elimination algorithm in modern one-based Fortran form, with native permutation, sparse matrix, regularization, and iterative-refinement wrappers.

The LDL package is LGPL-2.1-or-later.  Attribution and a copy of LGPL-2.1 are retained in `NOTICE.md` and `LICENSES/`.

## MatrixExtra integration

The optional adapter vendors the previously generated `MatrixExtra-fortran-v0.1.0` package under `integration/matrixextra-adapter/vendor/`.  MatrixExtra is used only for interoperability/construction, not for KKT factorization.

## Translation provenance

The Fortran source in `src/` is a language translation/reimplementation of the computational model and solver algorithms.  It does not embed the original ECOS C implementation bodies.

## v0.4 ordering/equilibration/certificate work

The v0.4 approximate-minimum-degree ordering, cone-preserving equilibration, sparse homogeneous ray-certificate builders, persistent workspace cache, and sparse ECOS_BB node-reuse code are new Fortran implementations in this translation. They are algorithmically motivated by the behavior required for ECOS-style sparse conic solving but are not copied SuiteSparse AMD or ECOS C source bodies.
