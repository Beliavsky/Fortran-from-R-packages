# ECOS / MatrixExtra adapter

Optional zero-densification interoperability layer between `ECOSolveR-fortran` and the previously translated `MatrixExtra-fortran` / `Matrix-fortran` sparse types.

The core ECOS solver does **not** depend on MatrixExtra. This project is useful when an application already represents matrices using `matrix_sparse:csc_matrix` or constructs them through MatrixExtra COO/CSR utilities.

Provided routines:

```text
ecos_csc_from_matrix
ecos_csc_from_coo
setup_problem_matrixextra
```

The adapter copies sparse indices/values directly into ECOS CSC storage. It deliberately does not route KKT factorization through MatrixExtra because the ECOS v0.4 core has its own sparse symbolic/numeric LDL backend, ordering, refinement, and cross-solve cache.

Build with FPM from this directory:

```text
fpm build
fpm run --example matrixextra_sparse_lp
```

or use `scripts/test_gfortran.sh` for the strict standalone validation build.
