# quantreg-fortran

Modern Fortran/FPM translation of computational portions of the R package
**quantreg 6.1**.

The public API works directly with numeric arrays. R formula parsing, model
frames, S3 methods, printing, graphics and dynamic-library registration are
intentionally not reproduced.

## Implemented numerical functionality

- Frisch-Newton interior-point quantile regression (`rq_fit_fnb`).
- Linear inequality constrained QR, with `R*b >= r` (`rq_fit_fnc`).
- Multiple-quantile Frisch-Newton fits (`rq_fit_qfnb`).
- Weighted QR (`rq_wfit_fnb`).
- Lasso-penalized quantile regression (`rq_fit_lasso`).
- Iteratively reweighted SCAD quantile regression (`rq_fit_scad`).
- Portnoy-Koenker preprocessing for large dense problems (`rq_fit_pfn`).
- Local linear quantile regression (`lprq`).
- Nonlinear QR by iterative linearization (`nlrq_fit`).
- Hyndman-Fan quantiles types 1--9 (`kuantiles`) and selection (`qselect`).
- Quantile check loss utilities.
- XY bootstrap for QR coefficients.
- Recursive least squares.
- Combination generation and exponential RNG utility.

The dense Frisch-Newton predictor/corrector algorithm is a direct translation
of the package's `rqfnb.f` / `rqfnc.f` logic. The translated default path is
self-contained and uses a small modern Cholesky solver instead of external
BLAS/LAPACK calls.

## Build

```text
fpm build
fpm test
```

Strict GNU Fortran validation scripts are also supplied in `scripts/`.

## Example

```fortran
use quantreg, only : dp, rq_result, rq_fit_fnb
real(dp) :: x(5,2), y(5)
type(rq_result) :: fit

x(:,1) = 1.0_dp
x(:,2) = [0.0_dp, 1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp]
y = [1.0_dp, 2.0_dp, 4.0_dp, 5.0_dp, 7.0_dp]
call rq_fit_fnb(x, y, 0.5_dp, fit)
```

See `API.md` and `TRANSLATION_COVERAGE.md` for details and explicitly listed
remaining gaps.
