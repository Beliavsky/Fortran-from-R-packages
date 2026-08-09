# API coverage map

This map relates major Matrix package concepts to the Fortran interfaces. R class dispatch is replaced by explicit array or derived-type procedures.

## Constructors and coercion

| Matrix concept | Fortran interface | Coverage |
|---|---|---|
| `Matrix`, dense coercion | ordinary `real(dp) :: a(:,:)` | Native Fortran representation |
| `sparseMatrix`, `spMatrix`, triplets | `csr_from_triplet` | Direct computational equivalent; duplicates summed |
| dense/sparse coercion | `csr_from_dense`, `csr_to_dense`, `csc_from_csr`, `csr_from_csc` | Direct |
| `Diagonal`, `sparseDiagonal` | `diag_matrix`, `sparse_diagonal`, `sparse_identity` | Direct numeric form |
| `band`, `bandSparse` | `band_matrix`, `csr_band` | Direct |
| `bdiag` | `block_diag` | Equal-sized dense blocks in a rank-3 array |
| `Hilbert` | `hilbert_matrix` | Direct |
| `rsparsematrix` | `random_sparse_matrix` | Direct numeric form |
| permutation matrices | `permutation_matrix`, `dense_permute`, `csr_permute` | Direct |
| packed/unpacked triangular | `pack_triangular`, `unpack_triangular` | Direct |

## Arithmetic and summaries

| Matrix concept | Fortran interface | Coverage |
|---|---|---|
| `%*%` | intrinsic `matmul`, `csr_matvec`, `csr_matmat`, `csr_multiply` | Direct |
| `crossprod`, `tcrossprod` | `crossprod`, `tcrossprod`, `csr_crossprod`, `csr_tcrossprod` | Direct |
| `kronecker` | `kronecker_product`, `csr_kronecker` | Direct |
| `KhatriRao` | `khatri_rao` | Direct |
| `rowSums`, `colSums`, means | dense and CSR summary functions | Direct |
| `nnzero`, `drop0` | `nnzero_dense`, `csr_matrix%nnz`, `csr_drop0` | Direct |
| `symmpart`, `skewpart` | same names | Direct |
| `isSymmetric`, triangular/diagonal tests | `is_symmetric`, `csr_is_symmetric`, `is_triangular`, `is_diagonal` | Direct numeric tests |
| `norm` | dense and CSR norm procedures | One, infinity, Frobenius, max-absolute |

## Factorizations and solves

| Matrix concept | Fortran interface | Coverage |
|---|---|---|
| `lu`, `solve`, `determinant` | `lu_factor`, `lu_solve`, `solve_linear`, `inverse_matrix`, `determinant`, `log_determinant` | Direct dense implementation |
| `chol`, positive-definite solve | `cholesky_factor`, `cholesky_solve` | Direct dense implementation |
| `qr`, least squares | `qr_factor`, `least_squares` | Direct modified Gram-Schmidt implementation |
| `BunchKaufman` | `ldlt_factor`, `ldlt_solve` | Unpivoted LDLT subset; not Bunch-Kaufman compatible |
| sparse LU/Cholesky/QR solves | `csr_lu_*`, `csr_cholesky_*`, `csr_least_squares` | Correct dense-backed adapters, not sparse symbolic algorithms |
| `rankMatrix` | `rank_matrix` | SVD-based |
| `kappa`, `rcond`, `condest` | `condition_number_1`, `reciprocal_condition_1` | Exact inverse-based estimate for moderate dense matrices |
| Schur complement | `schur_complement` | Direct block formula |
| general Schur decomposition | none | Not implemented |

## Matrix functions

| Matrix concept | Fortran interface | Coverage |
|---|---|---|
| `expm` | `matrix_exponential` | Scaling, Taylor series, and squaring |
| matrix powers | `matrix_power` | Integer powers, including negative powers |
| matrix square root | `matrix_sqrt_sym` | Symmetric positive-semidefinite matrices |
| `nearPD` | `near_positive_definite` | Alternating projections plus positive eigenvalue floor |
| positive-semidefinite projection | `projection_psd` | Direct |
| pseudoinverse | `pseudoinverse` | Compact-SVD based |
| null space | `null_space` | SVD/projector based |

## Sparse structure and I/O

| Matrix concept | Fortran interface | Coverage |
|---|---|---|
| compressed sparse classes | `csr_matrix`, `csc_matrix` | Numeric values, one-based canonical indices |
| `writeMM`, `readMM` | `write_matrix_market`, `read_matrix_market` | Coordinate real general/symmetric and pattern input |
| `readHB` | none | Not implemented |
| AMD/COLAMD | `minimum_degree_ordering`, `column_degree_ordering` | Portable greedy approximations, not SuiteSparse AMD/COLAMD |
| bandwidth reduction | `reverse_cuthill_mckee` | Direct graph algorithm |
| `dmperm` | none | Not implemented |

## R-only infrastructure omitted

S4 classes and validity methods, dimnames, coercion registries, model matrices, formulas, printing, plotting, localization, deprecated compatibility aliases, and external package adapters are outside the Fortran computational API.
