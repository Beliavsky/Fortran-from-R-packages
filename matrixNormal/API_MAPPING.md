# API mapping

This port translates the computational API of matrixNormal 0.1.2 into array-oriented modern Fortran.

| R function | Fortran API | Notes |
|---|---|---|
| `I(n)` | `identity_matrix(n)` | Renamed because Fortran is case-insensitive; a public `I` would collide with ordinary loop variable `i`. |
| `J(n,m)` | `ones_matrix(n,m)` | Renamed for the same reason (`J`/`j`). |
| `tr(A)` | `tr(A)` | Matrix trace. |
| `vec(A)` | `vec(A)` | Column-major vectorization, matching R's `as.vector(matrix)`. Names are R metadata and are omitted. |
| `vech(A)` | `vech(A)` | Correct lower-triangular half-vectorization in column-major order. |
| `is.square.matrix(A)` | `is_square_matrix(A)` | Logical scalar. |
| `is.symmetric.matrix(A,tol)` | `is_symmetric_matrix(A,tol)` | Uses the upstream sum-of-absolute-differences convention. |
| `is.positive.semi.definite(A,tol)` | `is_positive_semidefinite(A,tol)` | Symmetric-eigenvalue test with small values zeroed at `tol`. |
| `is.positive.definite(A,tol)` | `is_positive_definite(A,tol)` | Symmetric-eigenvalue test with upstream tolerance convention. |
| `dmatnorm(A,M,U,V,log)` | `dmatnorm(A,M,U,V,log_density,...)` | Uses Cholesky log determinants and SPD solves; no explicit inverses. |
| `pmatnorm(Lower,Upper,M,U,V,...)` | `pmatnorm(Lower,Upper,M,U,V,control,...)` | Returns `mvtnorm::probability_result` (value, error, status). Correct `V kron U` covariance is default. |
| `pmatnorm(-Inf,Upper,M,U,V,...)` | `pmatnorm(Upper,M,U,V,control,...)` | Four-argument generic overload. |
| `pmatnorm(-Inf,Inf,M,U,V,...)` | `pmatnorm(M,U,V,control,...)` | Three-argument full-support overload. |
| `rmatnorm(s=1,M,U,V,...)` | `rmatnorm(M,U,V,seed)` | One matrix draw. |
| `rmatnorm(s,M,U,V,...)` | `rmatnorm(s,M,U,V,seed)` | Extension: returns `n x p x s`; upstream currently ignores `s` and always draws one matrix. |

## Internal computational helpers

* Upstream `check_matnorm()` maps to public `check_matnorm()` with `ok` and `message` outputs instead of R exceptions.
* Upstream experimental `dmatnorm.logdet()` is not duplicated. The public Fortran `dmatnorm()` already uses stable Cholesky log determinants; the R helper is also affected by an expression-line break that prevents it from evaluating the full intended density formula.
* `find.eval()` is folded into the positive-(semi)definite predicates.
* R row/column names, warnings, printing, package-version checks, and data-frame conversion are interface behavior rather than numerical algorithms and are omitted.
