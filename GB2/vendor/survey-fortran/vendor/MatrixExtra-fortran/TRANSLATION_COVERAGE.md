# Translation coverage

## Translated computational families

The translation covers the main numerical behavior provided by MatrixExtra:

- CSR, CSC, COO and sparse-vector conversion paths.
- Sparse validation, index sorting, explicit-zero removal, filtering and mapping.
- CSR/COO arbitrary row/column slicing, including repeated selectors.
- Sparse scalar/row/column/block assignment.
- CSR row/column binding.
- Dense-CSC, CSR-dense, CSR-vector, sparse-vector, sparse-sparse multiplication,
  crossproducts and tcrossproducts.
- Sparse addition/subtraction/Hadamard multiplication.
- Logical AND/OR on sparse patterns.
- R column-major vector recycling for sparse-preserving arithmetic.
- MatrixExtra's zero-preserving mathematical functions.
- CSR norms and diagonal get/set.
- Opposite-format transpose conversion.

## Representation differences

R MatrixExtra specializes many S4 classes (`dgRMatrix`, `lgRMatrix`,
`ngRMatrix`, `TsparseMatrix`, triangular/symmetric classes, sparseVector,
float32, etc.). The Fortran API is numerical and array-oriented:

- numeric CSR/CSC use Matrix-fortran's `csr_matrix` / `csc_matrix`;
- COO and sparse vectors are added here;
- logical/pattern objects are represented as 0/1-valued sparse matrices;
- separate integer/logical/float32 class hierarchies are not reproduced.

## Performance differences

MatrixExtra's R `t_shallow()` can transpose CSR<->CSC by changing S4 metadata
without copying the data slots. The existing Matrix-fortran types own allocatable
arrays, so `csr_transpose_shallow` / `csc_transpose_shallow` copy the pointer,
index, and value arrays while preserving the same storage interpretation. A
future pointer/view type could make this genuinely zero-copy.

Core slicing, assignment, binding, elementwise arithmetic, and multiplication
remain sparse. Some highly specialized C++ kernels and OpenMP paths are not
ported one-for-one, so this release should not be expected to match MatrixExtra's
multithreaded performance.

## R-runtime behavior intentionally not reproduced

- S4 method registration and mutation of Matrix's global method behavior.
- `show`/`print` formatting and package options.
- R `NA` index semantics and all implicit-zero interactions with NA/NaN/Inf.
- `float::float32` overloads.
- dimnames, triangular/symmetric S4 metadata, and unit-diagonal slots.
- arbitrary repeated-index assignment's exact sequential R semantics; assignment
  selectors should be unique in the Fortran API.

The original sources are preserved under `original/MatrixExtra-master/` for
future parity work.
