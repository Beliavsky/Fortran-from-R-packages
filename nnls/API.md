# API

The convenience module is `nnls`.

## Kinds and status codes

```fortran
integer, parameter :: dp
integer, parameter :: NNLS_SUCCESS = 1
integer, parameter :: NNLS_BAD_DIMENSIONS = 2
integer, parameter :: NNLS_ITERATION_LIMIT = 3
```

## Result type

```fortran
type(nnls_result)
   real(dp), allocatable :: x(:)
   real(dp), allocatable :: fitted(:)
   real(dp), allocatable :: residuals(:)
   real(dp), allocatable :: dual(:)
   integer, allocatable :: passive(:)
   integer, allocatable :: bound(:)
   real(dp) :: rnorm
   real(dp) :: deviance
   integer :: mode
   integer :: nsetp
   integer :: iterations
end type
```

`passive` contains coefficient indices not bound at zero. `bound` contains the
remaining indices. `nsetp = size(passive)`.

## NNLS

```fortran
call nnls_fit(a, b, result [, max_iter])
```

`a` has shape `(m,n)` and `b` has length `m`. The returned solution satisfies
`x >= 0`.

## Mixed nonnegative/nonpositive least squares

```fortran
call nnnpls_fit(a, b, con, result [, max_iter])
```

`con` has length `n`. If `con(j) < 0`, coefficient `j` is constrained
nonpositive; otherwise it is constrained nonnegative.

## Generic name

The generic `nnls_solve` resolves to `nnls_fit` or `nnnpls_fit` according to
whether the `con` argument is supplied through the matching procedure form.
