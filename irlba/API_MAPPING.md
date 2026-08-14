# API mapping

| R `irlba` API | Fortran API | Notes |
|---|---|---|
| `irlba(A, nv, nu, ...)` | `irlba_svd(A, nv, nu, ...)` | Generic for dense real and `csc_matrix`. |
| restart with `v=<old result>` | `restart_from=old_result` | Reuses old U/V/d as a thick restart. |
| `scale=` | `scale=` | Implicit column scaling in matrix products. |
| `center=` | `center=` | Implicit column centering; preserves sparse products. |
| `shift=` | `shift=` | Scalar square-matrix shift. |
| `smallest=TRUE` | `smallest=.true.` | Full-SVD fallback in v0.1.0. |
| complex `irlba` | `irlba_complex` | Full complex LAPACK SVD compatibility path. |
| custom matrix class / `%*%` | extend `linear_operator`, call `irlb_operator` | Native matrix-free hook. |
| `partial_eigen` | `partial_eigen` | Dense and CSC overloads. |
| `prcomp_irlba` | `prcomp_irlba` | Dense truncated PCA. |
| `summary.irlba_prcomp` | fields `proportion`, `cumulative`, `total_variance` | Stored directly in `pca_result`. |
| `ssvd` | `ssvd` | Dense Shen-Huang sparse SVD/PCA. |
| `svdr` | `svdr` | Dense and native CSC sparse randomized SVD. |

R-specific S3/S4 dispatch, `Matrix` class conversion, warnings/messages,
`fastpath` selection, and deprecated R callback syntax are not translated as
runtime infrastructure.
