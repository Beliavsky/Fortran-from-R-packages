# CLA-fortran

A modern Fortran/FPM translation of the computational algorithms in R package
**CLA 0.96-3**, which implements Markowitz's Critical Line Algorithm for
box-constrained mean-variance portfolio optimization.

## Implemented numerical functionality

- Complete Critical Line Algorithm turning-point recursion
- Individual lower and upper asset-weight bounds
- Portfolio means and standard deviations for every turning point
- Frontier interpolation by target expected return (`find_sigma` / `findSig`)
- Frontier inversion by target standard deviation (`find_mu` / `findMu`)
- GARCH-derived expected-return and covariance inputs (`mu_sigma_garch`)
- Arbitrary GARCH(p,q) orders with normal or standardized Student-t innovations
- Original-name compatibility module `cla_api`
- Typed results, explicit status codes, and free-set masks

The R plotting and S3-printing methods are presentation infrastructure and are
not reproduced. Their numerical inputs are available directly in
`cla_result_t`.

## Build with FPM

```sh
fpm build
fpm test
fpm run
fpm run --example frontier_query
fpm run --example garch_inputs
```

BLAS and LAPACK are required. With GNU Fortran on systems without FPM:

```sh
./scripts/test_gfortran.sh strict
./scripts/test_gfortran.sh optimized
```

On Windows with MinGW-w64 or MSYS2, use the equivalent commands in
`scripts/test_gfortran.bat` after adjusting BLAS/LAPACK library names if needed.

## Basic use

```fortran
use kind_mod, only: dp
use cla, only: cla_result_t, critical_line

type(cla_result_t) :: result
result = critical_line(mu, covariance, lower_bounds, upper_bounds)
```

Original R-style names are available from `cla_api`:

```fortran
use cla_api, only: CLA, MS, findSig, findMu, muSigmaGarch
```

## Project layout

- `src/`: Fortran library
- `app/`: main demonstration
- `example/`: interpolation and GARCH examples
- `test/`: deterministic and structural tests
- `original/CLA/`: unchanged original package tree for provenance
- `PORTING.md`: compatibility and numerical differences
- `API.md`: public API inventory
- `TESTING.md`: validation details

## License

GPL-3.0-or-later, matching the original package metadata. See `LICENSE` and
`NOTICE.md`.
