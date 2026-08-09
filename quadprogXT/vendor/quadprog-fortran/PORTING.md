# Porting notes

## Preserved behavior

The numerical core remains the Goldfarb-Idnani dual active-set method used by
`quadprog` 1.5-8. The following upstream behavior is preserved:

- Objective convention `0.5*x^T*D*x - d^T*x`
- Constraint columns and equality-first ordering
- Factorized input convention (`R^{-1}` with `D=R^T R`)
- Lagrange multipliers and one-based active-set indices
- Major/deletion iteration counters
- Dense and native compact indexed-column algorithms
- Machine-precision feasibility safeguard added upstream for Talbot-Katz cases

## Modernization

- Fixed-form Fortran 77 was converted to free-form Fortran 2018.
- Labeled `do` loops were converted to structured named loops.
- LINPACK/BLAS calls were replaced with self-contained module procedures.
- The public API uses assumed-shape arrays, a derived result type, allocatable
  outputs, optional arguments, and explicit status constants.
- Inputs are copied because the upstream kernels overwrite their work arrays.
- Invalid inputs return status/message fields instead of raising R errors.

Some `goto` statements remain inside the private active-set state machine.
They preserve the control flow of the original, heavily tested algorithm and
are standard-conforming Fortran.

## Numerical assumptions

`D` must be symmetric positive definite. As upstream, the factorization uses
the upper triangle. The wrapper verifies dimensions and finiteness but does not
symmetrize `D`.

For compact constraints, duplicate variable indices within one column are not
combined. Callers should provide each nonzero variable at most once, matching
the upstream representation.
