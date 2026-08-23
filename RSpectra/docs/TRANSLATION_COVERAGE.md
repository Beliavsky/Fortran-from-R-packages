# Translation coverage

| Upstream computational area | Fortran status | Notes |
|---|---|---|
| `eigs()` dense matrix | Implemented | ARPACK regular; LAPACK full/shift fallback |
| `eigs_sym()` dense matrix | Implemented | ARPACK Lanczos, real shift-invert |
| Sparse general/symmetric matrix | Implemented | Native CSR, unshifted iterative |
| Sparse shift-invert | Partial | CSR densified before shifted factorization/selection |
| Function/operator eigen interface | Implemented | `linear_operator` extension API |
| `which=LM/SM/LR/SR/LI/SI` | Implemented | General |
| `which=LM/SM/LA/SA/BE` | Implemented | Symmetric |
| Initial vector, `ncv`, `tol`, `maxitr`, `retvec` | Implemented | `eigs_opts` |
| Real symmetric shift | Implemented | Iterative ARPACK mode 3 |
| General real/complex shift | Implemented | Explicit matrices via dense LAPACK spectrum + transformed selection |
| `svds()` dense | Implemented | Tall/wide normal equations |
| `svds()` sparse/operator | Implemented | Requires both product and transpose-product |
| SVD center/scale vectors | Implemented | Operator and dense |
| SVD automatic center/scale | Implemented for dense | Column means/norms |
| R `Matrix` S3 classes/coercion | Omitted | Replaced by CSR type |
| R S3 dispatch/list formatting/warnings | Omitted | Non-numerical interface behavior |
| Automatic symmetry detection in generic `eigs()` | Omitted | Call `eigs_sym()` explicitly |
| Spectra C API ABI | Not reproduced | Native Fortran API instead |
