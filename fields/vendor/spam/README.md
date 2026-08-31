# spam-fortran

Modern free-form Fortran/FPM port of the computational code in R package **spam 2.11-4**.

The design is intentionally numerical rather than an imitation of R's S4 object system. Sparse matrices use a small CSR derived type (`csr_matrix`), Ng-Peyton factors use `spam_chol`, and the public `spam` module re-exports the computational API.

## Major functionality

- CSR construction from dense/triplet input, duplicate consolidation, dense conversion and validation
- transpose, addition/subtraction, scaling, Hadamard and sparse matrix multiplication
- matrix-vector/matrix-matrix products, Kronecker products, row/column binding and subsetting
- diagonal operations, row/column sums and means, norms, bandwidth, cleaning and symmetry tests
- finite differences, Toeplitz/circulant/block-diagonal matrices, crossproducts and triangular extraction
- block `gmult`, permutations, callback entry mapping/margin reduction, sample column covariance
- Euclidean, maximum, Minkowski and great-circle sparse nearest-distance matrices
- `cov.exp`, spherical/correlation, nugget, Wu 1/2/3, Wendland 1/2, Matern 1/2, 3/2, 5/2, general Matern and Furrer/INLA-range Matern
- RW1/RW2/RWn/seasonal, regular/irregular IGMRF and regular-lattice GMRF precision matrices
- exact Ng-Peyton sparse Cholesky with MMD/RCM/no ordering, numerical updates, forward/back/full solves, determinants and `chol2inv`
- bundled symmetric and nonsymmetric ARPACK eigensolvers
- multivariate Normal simulation from covariance, precision and canonical forms
- shifted/Kshirsagar multivariate t simulation
- linearly constrained Normal simulation and conditional Normal simulation
- Gaussian random field simulation at supplied or generated grid locations
- Gaussian negative twice log likelihood and bounded MLE for built-in spam covariance models or user-supplied covariance callbacks
- upstream SparseKit, Ng-Peyton/PCx, ARPACK, permutation and sparse numerical subroutines retained in modern free-form source

## Build

```text
fpm build
fpm test
fpm run --example basic
```

BLAS and LAPACK are required and are declared in `fpm.toml`.

## Deliberately omitted R-facing code

Plotting/display/grid routines, print/summary/head/tail presentation, S4/S3 classes and coercion to Matrix classes, R option/profiling helpers, package version helpers, Germany map presentation/data helpers, and MatrixMarket/Harwell-Boeing file-interface wrappers are not translated. The native numerical algorithms used by the sparse matrix implementation are retained.

`apply.spam` is represented by typed callback entry mapping and scalar row/column reduction rather than R's dynamically typed arbitrary return objects. `spam64` is a separate R extension package and is not recreated here; the public CSR type currently uses the compiler's default integer kind, matching ordinary `spam` rather than `spam64`.

See `API_MAPPING.md`, `PORTING_NOTES.md`, `NOTICE.md`, and `VALIDATION.md`.
