# yieldcurves-fortran

A dependency-free modern Fortran/FPM translation of the numerical algorithms
in the R package **yieldcurves 0.1.0**.

The library fits and analyzes zero, par, and forward yield curves. It includes
Nelson-Siegel and Svensson models, cubic splines, interpolation, rate
conversions, fixed-income risk measures, carry and roll-down, and principal
component analysis.

## Build

```text
fpm build
fpm test
fpm run yieldcurves_demo
fpm run --example curve_fitting
fpm run --example risk_and_pca
```

The package uses only standard Fortran. No BLAS, LAPACK, optimization, or
statistics library is required.

## Main API

```fortran
use yieldcurves
```

Curve construction and fitting:

- `yc_curve`
- `yc_nelson_siegel`
- `yc_svensson`
- `yc_cubic_spline`
- `yc_fit`
- `yc_predict`
- `yc_interpolate`

Curve transformations and analysis:

- `yc_discount`
- `yc_forward`
- `yc_par_to_zero`
- `yc_zero_to_par`
- `yc_duration`
- `yc_bond_duration`
- `yc_zspread`
- `yc_key_rate_duration`
- `yc_carry`
- `yc_slope`
- `yc_level_slope_curvature`
- `yc_pca`

## Example

```fortran
program example
  use yieldcurves
  implicit none

  real(dp), parameter :: maturities(7) = [ &
    0.25_dp, 0.5_dp, 1.0_dp, 2.0_dp, 5.0_dp, 10.0_dp, 30.0_dp ]
  real(dp), parameter :: rates(7) = [ &
    0.052_dp, 0.050_dp, 0.048_dp, 0.045_dp, 0.042_dp, 0.040_dp, 0.043_dp ]

  type(curve_t) :: fit
  type(series_t) :: prediction

  fit = yc_nelson_siegel(maturities, rates)
  if (.not. fit%ok) error stop trim(fit%message)

  prediction = yc_predict(fit, [3.0_dp, 7.0_dp, 15.0_dp])
  print *, prediction%y
end program example
```

All public result structures contain `ok` and `message` fields. This replaces
R exceptions with explicit status handling suitable for libraries and batch
programs.

## Numerical implementation

- Nelson-Siegel fitting profiles out the three linear beta coefficients and
  performs bounded one-dimensional optimization over tau.
- Svensson fitting profiles out four linear beta coefficients and optimizes
  the two decay constants with bounded multi-start Nelder-Mead.
- Natural and FMM splines are implemented directly. The FMM boundary
  conditions are the not-a-knot conditions used by R's `splinefun`.
- PCA uses a symmetric Jacobi eigensolver on the covariance or correlation
  matrix and returns loadings, scores, standard deviations, and explained
  variance.
- Z-spreads use bracket expansion and bisection.
- The package uses double precision throughout.

## Scope

Plotting, S3 printing, CLI formatting, dates, data frames, and other R-specific
presentation infrastructure are not compiled. Every exported numerical
routine is represented by a Fortran procedure. See `COVERAGE.md` for the exact
mapping.

## License

MIT, preserving the license and copyright of the original package. The
unmodified upstream source is retained in `original/yieldcurves-0.1.0`, and
the supplied archive is retained in `provenance`.
