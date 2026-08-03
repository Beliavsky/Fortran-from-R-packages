# HierPortfolios-fortran

A modern Fortran translation of the computational code in the R package
`HierPortfolios` 1.0.2. The project is organized for the Fortran Package
Manager (FPM) and has no external numerical dependencies.

## Implemented strategies

- `HRP_Portfolio`: hierarchical risk parity.
- `HCAA_Portfolio`: hierarchical clustering-based asset allocation.
- `HERC_Portfolio`: hierarchical equal risk contribution.
- `DHRP_Portfolio`: constrained divisive hierarchical risk parity.

All routines accept a covariance matrix and return a `portfolio_result` with
weights, hierarchy order, cluster labels where applicable, status information,
and gap-statistic diagnostics when automatic cluster selection is used.

## Build

```sh
fpm build
fpm test
fpm run
```

GNU Fortran validation without FPM is available in `scripts/`.

## Minimal example

```fortran
use hierportfolios, only: dp, portfolio_result, HRP_Portfolio
real(dp) :: covar(4, 4)
type(portfolio_result) :: result

! Fill covar with a finite symmetric covariance matrix.
call HRP_Portfolio(covar, result, linkage='single')
if (result%ok()) print *, result%weights
```

## Scope

The numerical portfolio-allocation strategies are translated. R plotting,
S3/data-frame presentation, and bundled `.RData` loading are not part of the
Fortran library. The original package sources and data are retained under
`original/` for provenance.

See `PORTING.md` for numerical adaptations and `API.md` for procedure details.
