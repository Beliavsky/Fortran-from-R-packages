# MatrixExtra-fortran

Modern Fortran/FPM translation of the computational sparse-matrix functionality
of R package **MatrixExtra 0.1.15**.

The package builds on the previously translated `Matrix-fortran` CSR/CSC types
and adds MatrixExtra-style COO/sparse-vector handling, slicing and assignment,
binding, arithmetic, matrix multiplication, norms, diagonal operations, fast
storage-order transpose conversion, filtering/mapping, and R-style vector
recycling over sparse entries.

## Build

```text
fpm build
fpm test
```

No R, Rcpp, BLAS, or OpenMP dependency is required by the translated layer.
The vendored Matrix-fortran dependency is self-contained.

## Main types

```fortran
use matrixextra

type(csr_matrix)    :: a
type(csc_matrix)    :: ac
type(coo_matrix)    :: at
type(sparse_vector) :: v
```

`csr_matrix` and `csc_matrix` come from Matrix-fortran. COO matrices and sparse
vectors are added by MatrixExtra-fortran.

Logical and pattern sparse matrices are represented as ordinary sparse matrices
whose stored values are 0/1. `csr_patternize`, `csr_logical_and`, and
`csr_logical_or` provide pattern operations.

See `API.md` and `TRANSLATION_COVERAGE.md` for details.
