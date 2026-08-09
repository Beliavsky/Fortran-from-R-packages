# API

## Module `dykstra`

### Kind

```fortran
integer, parameter :: dp = kind(1.0d0)
```

### Result type

```fortran
type :: dykstra_result
    real(dp), allocatable :: solution(:)
    real(dp), allocatable :: unconstrained(:)
    real(dp) :: value
    integer :: iterations
    logical :: converged
    integer :: status
    character(len=:), allocatable :: message
end type
```

`value` is NaN when `factorized=.true.`, matching the R package's use of
`NA` in that mode.

### Solver

```fortran
call dykstra_solve(dmat, dvec, amat, result, &
    bvec=bvec, meq=meq, factorized=factorized, maxit=maxit, eps=eps)
```

Required arguments:

- `dmat(:,:)`: quadratic matrix `D`, or `R^{-1}` when factorized;
- `dvec(:)`: linear coefficient `d`;
- `amat(:,:)`: constraint matrix whose columns are constraint normals;
- `result`: `dykstra_result` output.

Optional arguments:

- `bvec(:)`: right-hand side; defaults to zero;
- `meq`: number of leading equality constraints; defaults to zero;
- `factorized`: defaults to false;
- `maxit`: maximum projection cycles; defaults to `30*n`;
- `eps`: convergence tolerance; defaults to `n*epsilon(1.0_dp)`.

Constraint convention:

```text
transpose(amat) * x >= bvec
```

The first `meq` rows of that system are interpreted as equalities.

### Status values

- `0`: converged;
- `1`: maximum number of cycles reached;
- negative values: invalid dimensions/options, non-PSD input, failed
  regularization, or a numerically zero constraint column.
