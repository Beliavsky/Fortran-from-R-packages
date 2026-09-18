# pbs

Modern free-form Fortran translation of the computational functionality in the R package **pbs** 1.1 by Shuangcai Wang.

The upstream package constructs periodic B-spline bases and provides an S3 prediction method that reuses the fitted spline specification. This translation preserves those numerical operations in a small dependency-free FPM package.

## Features

- Periodic polynomial B-spline bases with folded endpoint basis functions.
- Ordinary nonperiodic B-spline bases.
- R-compatible type-7 quantile selection of internal knots when `df` is supplied.
- Optional intercept column.
- Explicit boundary knots or boundaries inferred from the non-missing predictor range.
- IEEE NaN rows preserved as missing output rows.
- R-style nonperiodic extrapolation using a boundary Taylor expansion through the spline degree.
- Stored spline specification and `predict_pbs` evaluation on new predictor values.
- No BLAS, LAPACK, or external package dependency.

## Basic use

```fortran
use pbs_mod, only : dp, pbs_basis, pbs

type(pbs_basis) :: b
real(dp) :: x(5)

x = [0.0_dp, 0.5_dp, 1.0_dp, 1.5_dp, 2.0_dp]
b = pbs(x, knots=[0.5_dp, 1.0_dp, 1.5_dp], degree=3, &
   boundary_knots=[0.0_dp, 2.0_dp], periodic=.true.)
```

The returned `pbs_basis` object contains `values`, the internal knots, boundary knots, degree, intercept flag, periodic flag, and a status/message pair for argument errors.

## Translation coverage

Package status: **substantial**.

**2 of 2 (100%)** exported computational R functions/methods are mapped. This fraction measures mapped computational functions, not complete compatibility with the R object system.

| R function | Fortran module | Fortran procedure(s) | Mapping status |
|---|---|---|---|
| `pbs` | `pbs_api` | `pbs` | substantial |
| `predict.pbs` | `pbs_api` | `predict_pbs` | substantial |

`makepredictcall.pbs` is not counted because it rewrites R formula calls and does not perform a numerical spline operation. Plotting, formula processing, R class construction, dimnames, and S3 dispatch are intentionally omitted.

The numerical spline construction is independently checked against SciPy's `BSpline` implementation. On the permanent parity cases the maximum absolute differences are approximately `1.1e-16` for the periodic basis and `4.0e-15` for the nonperiodic extrapolation case.

## Build

With FPM installed:

```text
fpm build
fpm test
```

The package is written entirely in free-form Fortran and uses one real kind, `dp`, defined in `pbs_kinds` from `iso_fortran_env::real64`.

## Upstream and license

The upstream `DESCRIPTION`, `NAMESPACE`, R implementation, and manual page are retained under `upstream/` for provenance. The upstream package is GPL-2; see `LICENSE` and `NOTICE.md`.
