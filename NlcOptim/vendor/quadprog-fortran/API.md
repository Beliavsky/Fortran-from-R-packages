# API

## Modules

### `quadprog_kinds`

```fortran
integer, parameter :: dp = kind(1.0d0)
```

### `quadprog`

Exports `qp_result`, `solve_qp`, `solve_qp_compact`, and status constants.

## Result type

```fortran
type :: qp_result
  real(dp), allocatable :: solution(:)
  real(dp), allocatable :: unconstrained_solution(:)
  real(dp), allocatable :: lagrangian(:)
  integer, allocatable :: active_set(:)
  real(dp) :: value
  integer :: iterations(2)
  integer :: status
  character(len=:), allocatable :: message
contains
  procedure :: succeeded
end type
```

`iterations(1)` is the number of major active-set iterations.
`iterations(2)` is the number of previously active constraints deleted.
Constraint indices are one-based, as in R and Fortran.

## `solve_qp`

```fortran
fit = solve_qp(dmat, dvec, amat [, bvec] [, meq] [, factorized])
```

Arguments:

- `dmat(n,n)`: positive-definite quadratic matrix `D`, or `R^{-1}` when
  `factorized=.true.` and `D=R^T R`.
- `dvec(n)`: linear coefficient `d`.
- `amat(n,q)`: constraint normals as columns.
- `bvec(q)`: right-hand sides; defaults to zero.
- `meq`: number of leading equality constraints; defaults to zero.
- `factorized`: whether `dmat` contains upper-triangular `R^{-1}`.

The objective is `0.5*x^T*D*x - d^T*x`. The first `meq` constraints are
equalities; remaining constraints are `amat(:,j)^T*x >= bvec(j)`.

## `solve_qp_compact`

```fortran
fit = solve_qp_compact(dmat, dvec, values, indices &
  [, bvec] [, meq] [, factorized])
```

For `q` constraints and at most `m` nonzeros per constraint:

- `values(m,q)` stores nonzero coefficients.
- `indices(m+1,q)` stores the count in row 1 and variable indices in rows
  `2:count+1`.

For constraint `j`, coefficient `values(k,j)` belongs to variable
`indices(k+1,j)`. This is the same compact representation used by the R
package. The solver works directly on this representation and does not form a
dense `n x q` constraint matrix.

## Status constants

- `qp_success = 0`
- `qp_inconsistent_constraints = 1`
- `qp_not_positive_definite = 2`
- `qp_invalid_dimensions = 3`
- `qp_invalid_meq = 4`
- `qp_invalid_compact_index = 5`
- `qp_nonfinite_input = 6`

Call `fit%succeeded()` to test for `qp_success`.
