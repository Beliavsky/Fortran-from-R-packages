# Translation notes

## Upstream

- Package: `tensorA`
- Version: 0.36.2.1
- Author: K. Gerald van den Boogaart
- License: GPL (>= 2)
- Upstream compiled code: `src/tensora.c`

The upstream C file contains two batched matrix-multiplication helpers, one
real and one complex. The remainder of the package is R code. The Fortran
translation does not compile or call that C code.

## R-to-Fortran mapping

| Upstream operation | Fortran API |
|---|---|
| `to.tensor`, `as.tensor` | `tensor`, `to_tensor`, `tensor_from_real`, `tensor_from_complex` |
| `norm.tensor` | `norm_tensor`, `norm_tensor_along` |
| `opnorm.tensor` | `opnorm_tensor`, `opnorm_by_tensor` |
| `margin.tensor` | `margin_tensor` |
| `diagmul.tensor` | `diagmul_tensor` |
| `pos.tensor` | `pos_tensor` |
| `reorder.tensor` | `reorder_tensor`, `reorder_tensor_pos`, `reorder_tensor_names` |
| `mul.tensor` | `mul_tensor`, `mul_tensor_pos`, `mul_tensor_names` |
| `rep.tensor` | `repeat_tensor` |
| `trace.tensor` | `trace_tensor` |
| `delta.tensor` | `delta_tensor` |
| `diag.tensor` | `diag_tensor` |
| `tripledelta.tensor` | `tripledelta_tensor` |
| `one.tensor` | `one_tensor` |
| `mark.tensor` | `mark_tensor` |
| `inv.tensor` | `inv_tensor` |
| `solve.tensor` | `solve_tensor` |
| `chol.tensor` | `chol_tensor` |
| `svd.tensor` | `svd_tensor` |
| `power.tensor` | `power_tensor` |
| `to.matrix.tensor` | `to_matrix_tensor` |
| `untensor` | `untensor_tensor` |
| `slice.tensor` | `slice_tensor` |
| `undrop.tensor` | `undrop_tensor` |
| `bind.tensor` | `bind_tensor` |
| `toPos.tensor` | `positions_by_name`, `names_to_positions` |
| `einstein.tensor`, `%e%` | `einstein_pair` |
| tensor `+ - * /` | overloaded operators and `add/sub/elem_mul/elem_div_tensor` |
| `contraname` | `contraname` |
| `is.covariate`, `is.contravariate` | `is_covariate_name/tensor`, `is_contravariate_name/tensor` |
| `as.covariate`, `as.contravariate` | `as_covariate_name`, `as_contravariate_name` |
| `drag.tensor` | `drag_tensor` |
| `riemann.tensor`, `%r%` | `riemann_pair` |
| `mean.tensor` | `mean_tensor` |
| `var.tensor` | `var_tensor`, `cov_tensor` |

## Numerical implementation

### Contraction

The original C backend treats each contraction as a sequence of matrix
multiplications over optional parallel dimensions. The Fortran implementation
uses the same organization: reorder axes, reshape into matrix batches, then
use intrinsic `matmul` for each batch. Both real and complex tensors use this
native path.

### Linear algebra

The package is self-contained. Complex one-sided Jacobi SVD, Gauss-Jordan
inverse, SVD-based Moore-Penrose pseudoinverse and Cholesky factorization are
implemented in Fortran. This avoids a mandatory BLAS/LAPACK dependency while
retaining complex tensor support.

### Storage

Fortran column-major storage matches R array storage. This makes reshape-based
matrix blocks natural and avoids changing the original tensor element order.

## Intentional interface differences

The computational tensor algebra is translated, but R-specific object syntax
is not reproduced literally:

1. R S3 methods (`$`, `^`, `|`, replacement methods, `ftable`) are replaced by
   explicit procedures such as `rename_axis`, `rename_first_axes`,
   `reorder_tensor`, and `slice_tensor`.
2. `tensor_t` stores named axes but not R `dimnames` level labels. Slicing is by
   integer position rather than character level label.
3. Variadic `einstein.tensor(...)`/`riemann.tensor(...)` calls with embedded R
   rename/diagonal instructions are expressed as pairwise `einstein_pair` or
   `riemann_pair` calls plus explicit `rename_axis`, `mark_tensor` or
   `diagmul_tensor` operations.
4. R's flexible list-driven `untensor` can be composed from repeated
   `untensor_tensor` calls in Fortran.
5. `var.tensor(..., na.rm=TRUE)` pairwise-complete missing-value semantics are
   not reproduced. `mean_tensor` supports optional finite-value removal;
   covariance expects finite data.
6. Printing/formatting and R class metadata are intentionally omitted.

## Source form and dependencies

All compiled Fortran sources are free-format `.f90`. The FPM library does not
compile or link the original C source retained in `upstream/`.
