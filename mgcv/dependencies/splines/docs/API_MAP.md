# R to Fortran API map

| Original R/C routine | Fortran API | Notes |
|---|---|---|
| `spline.des`, `splineDesign`, `spline_basis` | `spline_design` | Dense design matrix; per-point derivative orders supported. |
| `spline.value`, `predict.bSpline`, `spline_value` | `b_spline_t%evaluate` | Scalar and vector evaluation; derivative order optional. |
| `bs` | `bs_basis` | Explicit knots or R type-7 quantile knots selected from `df`; polynomial extrapolation outside boundaries. |
| `ns` | `natural_spline_basis` | Natural cubic basis with linear extrapolation and two endpoint constraints. |
| `interp.spline`, `interpSpline.default` | `fit_interpolating_spline` | Natural cubic interpolation returned as `b_spline_t`. |
| `periodic.spline`, `periodicSpline.default` | `fit_periodic_spline` | Arbitrary spline order and optional custom knot sequence. |
| `polySpline.bSpline` | `to_polynomial_spline` | Converts to local Taylor coefficients in `poly_spline_t`. |
| `predict.polySpline` | `poly_spline_t%evaluate` | Horner evaluation and derivatives. |
| `predict.nbSpline`, `predict.npolySpline` | Natural flags on typed objects | Linear extrapolation is performed by the object evaluator. |
| `predict.pbSpline`, `predict.ppolySpline` | Periodic flags on typed objects | Inputs are wrapped into one period before evaluation. |
| `backSpline` | `inverse_monotone_spline` | Inverts a monotone cubic polynomial spline. |
| `linear.interp`, `lin_interp` | `linear_interp` | Sorts input pairs and supports arbitrary query order. |
| `splineKnots`, `splineOrder`, `coef` | Public type components | `knots`, `order`, and `coefficients`. |
| `predict.bs`, `predict.ns` | Re-run the corresponding basis routine | Supply the retained knots, boundaries, degree, and intercept choice. |

## Intentionally omitted

- Formula methods and evaluation environments.
- S3 class dispatch and R object attributes.
- `xyVector` and data-frame conversion helpers.
- Plot and print methods.
- Dynamic-library registration and `.C` interface wrappers.

These omissions are interface or presentation code, not numerical algorithms.
