# desirability-fortran

A modern Fortran/FPM translation of the computational code in Max Kuhn's R package
`desirability` version 2.1.

The original R package is licensed GPL-2. This translation preserves that license;
see `COPYING`.

## Scope

Translated numerical behavior:

- `dMax` -> `d_max`
- `dMin` -> `d_min`
- `dTarget` -> `d_target`
- `dArb` -> `d_arb`
- `dBox` -> `d_box`
- `dCategorical` -> `d_categorical`
- `dOverall` -> `d_overall`
- all computational `predict.*` behavior
- non-informative missing-value desirabilities
- zero-desirability tolerance replacement
- arbitrary-function linear interpolation with endpoint carry-forward
- categorical desirabilities
- overall geometric-mean aggregation

Plotting code from `R/plot.R` is intentionally omitted. R-specific S3 dispatch, calls,
printing, data-frame labels, and graphics infrastructure are not translated.

## Build

```text
fpm build
fpm test
fpm run --example basic
```

The project requires Fortran 2008 or later. It has no external dependencies.

## Basic numeric use

```fortran
use desirability, only : dp, d_max, d_target, d_overall, &
    d_overall_type, hold, predict

type(d_overall_type) :: overall
real(dp) :: x(1, 2)
real(dp), allocatable :: score(:)

overall = d_overall([hold(d_max(80.0_dp, 97.0_dp)), &
    hold(d_target(55.0_dp, 57.5_dp, 60.0_dp))])
x(1, :) = [81.09_dp, 59.85_dp]
score = predict(overall, x)
```

`predict_all` returns the individual desirabilities plus the overall value in the last
column.

## Missing values

For numeric input, IEEE quiet NaN is treated like R's `NA`. By default it is replaced
by the object's non-informative desirability, computed with the same 100-point grid as
the R package. An explicit `missing=` argument can override this for standalone
predictions.

For categorical input, an empty/blank string represents a missing category because the
R package itself disallows empty category names.

## Mixed numeric/categorical overall desirability

Use `desirability_input` values made with `numeric_input`, `categorical_input`, and
`missing_input`, then pass the resulting rank-2 array to `predict` or `predict_all`.

## Notes on compatibility

The computational formulas and boundary conventions follow the R 2.1 source. `dArb`
uses linear interpolation between sorted knots and carries the extreme supplied
values beyond the input range. Tied interpolation x-values are averaged for the
interpolation itself, matching R's `approx` default tie handling; values beyond the
range retain the first/last sorted supplied desirability as in the package code.
