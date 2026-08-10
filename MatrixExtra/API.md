# API

## Conversions

- `coo_from_dense`, `coo_to_dense`
- `coo_from_csr`, `csr_from_coo`
- `coo_from_csc`, `csc_from_coo`
- `csr_to_csc`, `csc_to_csr`
- `sparse_vector_from_dense`, `sparse_vector_to_dense`
- `csr_from_sparse_vector`, `coo_from_sparse_vector`
- `csr_transpose_shallow`, `csc_transpose_shallow`

The last two preserve MatrixExtra's storage reinterpretation semantics (CSR of
`A` becomes CSC of `transpose(A)` with identical pointer/index/value arrays),
but Fortran allocatable assignment copies those arrays rather than sharing R
slots.

## Validation/utilities

- `csr_check`, `coo_check`, `sparse_vector_check`
- `csr_sort_indices`, `coo_sort_indices`, `sparse_vector_sort_indices`
- `csr_remove_zeros`, `coo_remove_zeros`, `sparse_vector_remove_zeros`
- `empty_sparse`
- `csr_filter`
- `csr_map`

## Slicing and assignment

- `csr_get`
- `csr_slice`
- `coo_slice`
- `csr_set_value`
- `csr_set_row_constant`
- `csr_set_col_constant`
- `csr_set_block_constant`
- `csr_assign_dense_block`

Slicing accepts arbitrary/repeated row and column selectors. Assignment expects
unique selectors; unlike R, NA indices are not a Fortran indexing concept.

## Binding

- `csr_rbind`
- `csr_cbind`

## Matrix/vector products

- `csr_matvec_extra`
- `csr_sparse_vector_matmul`
- `csr_dense_matmul`
- `dense_csc_matmul`
- `csr_csr_matmul`
- `csr_crossprod`
- `csr_tcrossprod`

## Sparse elementwise operations

- `csr_elem_add`
- `csr_elem_subtract`
- `csr_elem_multiply`
- `csr_logical_and`
- `csr_logical_or`
- `csr_patternize`
- `csr_scale_values`
- `csr_divide_scalar`
- `csr_power_scalar`

## R-style vector recycling

For a stored matrix entry `(i,j)`, the vector position is computed from the
column-major R index `i + (j-1)*nrow`, recycled modulo `size(v)`.

- `csr_multiply_vector`
- `csr_divide_vector`
- `csr_power_vector`
- `csr_mod_vector`
- `csr_intdiv_vector`

## Zero-preserving mathematical functions

- `csr_apply_sqrt`
- `csr_apply_abs`
- `csr_apply_log1p`
- `csr_apply_expm1`
- `csr_apply_sin`
- `csr_apply_sinh`
- `csr_apply_tan`
- `csr_apply_tanpi`
- `csr_apply_tanh`
- `csr_apply_atanh`
- `csr_apply_sign`
- `csr_apply_floor`
- `csr_apply_ceiling`
- `csr_apply_trunc`
- `csr_apply_round`
- `csr_apply_signif`

## Linear algebra

- `csr_norm` with norm types `1/O`, `I`, `F`, `M`, `2`
- `csr_diag`
- `csr_set_diag`

The spectral 2-norm uses power iteration rather than Matrix's external dense or
sparse singular-value machinery.
