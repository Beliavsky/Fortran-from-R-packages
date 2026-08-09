# API

## Modules

### `quadprogxt`

Primary translation module.

### `quadprog`

Provided by the vendored `quadprog-fortran` dependency.  Its public
`qp_result`, `solve_qp`, and `solve_qp_compact` types/routines are used directly.

## `solve_qp_xt`

```fortran
fit = solve_qp_xt(dmat, dvec [, amat] [, bvec] [, meq] [, factorized], &
  [, amat_posneg] [, bvec_posneg] [, dvec_posneg] [, b0], &
  [, amat_posneg_delta] [, bvec_posneg_delta] [, dvec_posneg_delta], &
  [, tol] [, compact] [, normalize])
```

The ordinary objective is

`0.5*b^T*Dmat*b - dvec^T*b`.

Constraint columns in `amat` mean

`amat(:,j)^T*b >= bvec(j)`,

with the first `meq` columns treated as equalities.

When absolute-value terms are requested, the internal decision vector is

`z = [b, abs_b_aux, abs_delta_aux]`.

The auxiliary constraints enforce

- `abs_b_aux >= abs(b)`, and
- `abs_delta_aux >= abs(b-b0)`.

`amat_posneg` acts on the conceptual vector
`[b_positive,b_negative]`, while `amat_posneg_delta` acts on
`[delta_positive,delta_negative]`.  The exact linear map used by the R package
is applied before solving.

`dvec_posneg` and `dvec_posneg_delta` have length `2*n` and provide the linear
objective coefficients in those conceptual positive/negative coordinates.

The return type is the vendored `quadprog::qp_result`.  The original decision
variables are always `fit%solution(1:n)`.

## `build_qp_xt`

Same arguments as `solve_qp_xt`; returns a `qpxt_problem` containing the
expanded ordinary quadratic program rather than solving it.

```fortran
type :: qpxt_problem
  real(dp), allocatable :: dmat(:,:)
  real(dp), allocatable :: dvec(:)
  real(dp), allocatable :: amat(:,:)
  real(dp), allocatable :: bvec(:)
  real(dp), allocatable :: compact_amat(:,:)
  integer, allocatable :: aind(:,:)
  integer :: meq
  logical :: factorized
  logical :: compact
  integer :: n_original
  integer :: n_abs
  integer :: n_delta
  integer :: status
  character(len=:), allocatable :: message
contains
  procedure :: succeeded
end type
```

## `normalize_constraints`

```fortran
normed = normalize_constraints(amat, bvec)
```

For each constraint column `j`, divides both the column and `bvec(j)` by
`sqrt(dot_product(amat(:,j),amat(:,j)))`.  A zero-norm constraint is rejected.

## `convert_to_compact`

```fortran
compact = convert_to_compact(amat)
```

For `q` constraints, returns:

- `compact%amat(max_nnz,q)`: nonzero coefficient values;
- `compact%aind(max_nnz+1,q)`: count in row 1, one-based variable indices in
  later rows.

This is the compact layout used by `quadprog::solve_qp_compact` and the R
package `quadprog::solve.QP.compact`.

## Translation status codes

- `qpxt_success = 0`
- `qpxt_invalid_dimensions = 100`
- `qpxt_zero_constraint = 101`
- `qpxt_missing_b0 = 102`
- `qpxt_nonfinite_input = 103`
- `qpxt_invalid_tolerance = 104`

Solver errors after a successful build use the status codes from the vendored
`quadprog` package.
