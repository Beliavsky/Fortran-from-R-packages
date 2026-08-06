# matrix-fortran

A modern Fortran 2018/FPM translation of the central computational ideas in the R package **Matrix 1.7-6**.

The port uses ordinary Fortran arrays for dense matrices and native derived types for compressed sparse row (CSR) and compressed sparse column (CSC) matrices. It is self-contained and does not require R, BLAS, LAPACK, CHOLMOD, or the other bundled SuiteSparse libraries.

## Implemented computational areas

### Dense matrices

- Identity, diagonal, banded, block-diagonal, Toeplitz, Hilbert, permutation, and companion constructors
- Kronecker and Khatri-Rao products
- Symmetric and skew-symmetric parts
- Row/column sums and means, scaling, crossproducts, packing, unpacking, and permutations
- Matrix norms, trace, structural tests, and nonzero counts
- LU factorization with pivoting, determinant, log-determinant, solves, and inverse
- Cholesky factorization and positive-definite solves
- QR factorization and least squares
- Symmetric eigendecomposition
- Compact singular-value decomposition, numerical rank, pseudoinverse, and null space
- LDLT factorization for nonsingular symmetric matrices
- Matrix exponential, integer matrix powers, symmetric square root, Schur complements
- Nearest positive-semidefinite/positive-definite and nearest-correlation calculations
- Exact one-norm condition number for moderate dense matrices

### Sparse matrices

- Canonical one-based CSR and CSC representations
- Construction from dense matrices or triplets, with duplicate aggregation and sorted indices
- Dense/CSR/CSC conversion and transpose
- Sparse matrix-vector and sparse-dense matrix products
- Sparse addition, multiplication, scaling, zero dropping, band extraction, diagonals, permutations, and Kronecker products
- Sparse row/column summaries, norms, trace, crossproducts, and symmetry checks
- Reverse Cuthill-McKee, greedy minimum-degree, and column-degree orderings
- Matrix Market coordinate read/write
- Dense-backed LU, Cholesky, and least-squares adapters for sparse inputs

## Important scope distinction

This is not a source-level replacement for the complete R package. Matrix contains a large S4 class hierarchy and approximately 165,000 lines of R, C, and bundled SuiteSparse code. The following are intentionally not reproduced:

- R S4 classes, method dispatch, coercion tables, dimnames, model-matrix integration, printing, localization, and plotting
- Native CHOLMOD/CXSparse/AMD/COLAMD/METIS bindings and their symbolic factor objects
- Supernodal sparse Cholesky updates/downdates and sparse factor reuse
- Exact Bunch-Kaufman, general real Schur, Dulmage-Mendelsohn, Harwell-Boeing, and sparse-QR compatibility interfaces
- Logical, pattern, integer, and complex class hierarchies as separate runtime types

Sparse factorizations in this release convert to dense storage. They are correct for moderate matrices but do not have SuiteSparse scalability. See `API_MAP.md` and `PORTING_NOTES.md`.

## Build with FPM

```text
fpm build
fpm test
fpm run demo_matrix
fpm run --example dense_solve
fpm run --example sparse_poisson
```

## Build without FPM

On Unix-like systems:

```text
./build_gfortran.sh
```

On Windows with `gfortran` on `PATH`:

```text
build_gfortran.bat
```

## Minimal example

```fortran
program example
   use matrix, only : dp, solve_linear, matrix_success
   implicit none
   real(dp), allocatable :: a(:,:), b(:,:), x(:,:)
   integer :: info

   a = reshape([3.0_dp, 1.0_dp, 1.0_dp, 2.0_dp], [2, 2])
   b = reshape([9.0_dp, 8.0_dp], [2, 1])
   call solve_linear(a, b, x, info)
   if (info /= matrix_success) error stop 'solve failed'
   print *, x(:, 1)
end program example
```

## Numerical design

- `dp` is defined as `kind(1.0d0)`.
- Sparse indices are one-based and each CSR row or CSC column is canonical and strictly sorted.
- Algorithms are written in portable Fortran rather than delegated to external libraries.
- The symmetric eigensolver uses Jacobi rotations; the compact SVD is constructed from the eigendecomposition of `transpose(A)*A`.
- `near_positive_definite` follows alternating projections with an eigenvalue floor.

For very large or ill-conditioned production problems, a BLAS/LAPACK/SuiteSparse backend would be substantially faster and often more accurate.

## License

GPL-3.0-only for the translated project. The complete upstream source and all bundled third-party license notices are retained under `original/`.
