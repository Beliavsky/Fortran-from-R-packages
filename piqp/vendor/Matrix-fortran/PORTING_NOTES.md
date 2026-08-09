# Porting notes

## Dense representation

Use ordinary rank-2 `real(dp)` arrays. Fortran is column-major, as is R, so dense storage order maps naturally.

## Sparse representation

`csr_matrix` stores:

- `nrow`, `ncol`
- `row_ptr(1:nrow+1)`
- `col_ind(1:nnz)`
- `values(1:nnz)`

Indices are one-based. Each row has strictly increasing column indices. `csr_from_triplet` sorts entries, sums duplicates, and removes values at or below the requested tolerance.

`csc_matrix` is analogous with `col_ptr` and `row_ind`.

## Error handling

Procedures return integer status values from `matrix_status`:

- `matrix_success`
- `matrix_err_shape`
- `matrix_err_singular`
- `matrix_err_not_posdef`
- `matrix_err_invalid`
- `matrix_err_convergence`
- `matrix_err_io`

## Numerical limitations

The implementation is intentionally self-contained. It does not call optimized BLAS/LAPACK or SuiteSparse routines. Consequences include:

- Dense algorithms are most suitable for small and medium matrices.
- Sparse factorization adapters densify their input.
- The Jacobi eigensolver is robust for symmetric matrices but slower than LAPACK.
- The SVD uses the normal matrix `A^T A`, which can square the condition number.
- The unpivoted LDLT routine is not a replacement for Bunch-Kaufman on difficult indefinite systems.
- Matrix exponential and nearest-positive-definite results can differ slightly from Matrix package results because the algorithms and stopping rules differ.

## Future extension points

The public API can be retained while replacing internals with:

- BLAS/LAPACK through standard interfaces
- a SuiteSparse or MUMPS backend for sparse factorization
- Householder QR and bidiagonal SVD
- pivoted symmetric-indefinite factorization
- complex and logical/pattern sparse value types
