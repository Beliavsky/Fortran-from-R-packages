# Translation coverage

## Translated

All exported computational functions are represented:

- `nnls`
- `pnnls`
- `pnnqp`
- `ldp`
- `lsi`
- `lsei`
- `hfti`
- `qp`
- `indx`
- `matMaxs`

The quadratic-programming transformations use a self-contained symmetric
Jacobi eigensolver.  Bounds are translated to inequality rows using the same
sign convention as the R package.

## Implementation differences

The modern code preserves the mathematical algorithms but not every mutated
work array of the 1970s Fortran interface.

- NNLS retains the Lawson-Hanson active/passive-set algorithm but recomputes a
  passive least-squares solve instead of updating Householder factors in place.
- Partial NNLS eliminates the unrestricted columns with a pseudoinverse/null
  projection before applying NNLS to the restricted columns.  This is
  mathematically equivalent but does not reproduce the original transformed
  `Q*A` work matrix byte-for-byte.
- `hfti_solve` returns the solution, rank and residual norms.  Its
  `transformed_b` field places the solution in the leading rows but does not
  attempt to reproduce the complete internal `Q*b` work array from HFTI.
- `lsi/lsei` use equality null-space reduction, LDP feasibility, and a feasible
  active-set quadratic step.  This is the same constrained least-squares
  problem as the R QR/SVD reductions but uses modern standalone linear algebra.

The complete original Lawson-Hanson/Yong Wang source is retained in
`original/lsei-master/` for exact provenance and further parity work.
