# API

All public definitions are available through:

```fortran
use hierportfolios
```

The real kind is `dp = kind(1.0d0)`.

## Result type

```fortran
type(portfolio_result)
  integer :: status
  character(len=:), allocatable :: message
  character(len=:), allocatable :: method
  integer :: n_clusters
  integer :: iterations
  real(dp), allocatable :: weights(:)
  integer, allocatable :: order(:)
  integer, allocatable :: clusters(:)
  real(dp), allocatable :: gap(:)
  real(dp), allocatable :: gap_se(:)
contains
  procedure :: ok
end type
```

Status constants are `hp_success`, `hp_invalid_argument`,
`hp_numerical_failure`, and `hp_not_converged`.

## HRP_Portfolio

```fortran
call HRP_Portfolio(covar, result [, linkage])
```

`linkage` may be `single`, `complete`, `average`, or `ward`. The default is
`single`. The routine converts covariance to correlation distance, clusters
the rows of that distance matrix, and applies recursive bisection using
inverse-variance cluster portfolios.

## HCAA_Portfolio

```fortran
call HCAA_Portfolio(covar, result [, linkage, clusters, gap_references, seed])
```

The default linkage is `ward`. Supply `clusters` to choose the terminal count
explicitly. Otherwise a deterministic gap-statistic simulation is used;
`gap_references` defaults to 100 and `seed` controls the portable generator.

## HERC_Portfolio

```fortran
call HERC_Portfolio(covar, result [, linkage, clusters, gap_references, seed])
```

The terminal clusters use normalized inverse-variance weights. Capital is
allocated through successive hierarchy splits in inverse proportion to the
risk of the two child groups.

## DHRP_Portfolio

```fortran
call DHRP_Portfolio(covar, result [, tau, ub, lb])
```

The routine uses DIANA divisive clustering. `tau` lies in `[0,1]` and controls
how far the split search may move from the midpoint. `lb` and `ub` contain
per-asset lower and upper bounds. Defaults are zero and one. Feasible bounds
must satisfy `sum(lb) <= 1 <= sum(ub)`.
