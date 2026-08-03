# Porting notes

## Object model

The upstream package uses R S4 classes that inherit from numeric and character
vectors. The translation uses explicit Fortran derived types. Operations that
R expresses through S3/S4 dispatch are exposed as typed procedures such as
`spot_compound`, `curve_discount`, and `forwardrate_from_curve`.

R dots are written as underscores in Fortran names. The replacement function
`interpolation<-` is represented by `set_interpolation`.

## Calendars

The upstream package delegates date arithmetic to `bizdays`, whose named
calendars contain external holiday data. This self-contained translation
supports:

- `actual`: calendar-day offsets;
- `weekdays` or `business`: Monday-Friday offsets without holidays.

Other names return `FI_UNSUPPORTED_CALENDAR`/`FI_INVALID_ARGUMENT`. They are
not silently treated as actual days. Day-count conversion itself is fully
implemented and does not require calendar data.

## Interpolation

Linear, log-linear, natural cubic, and flat-forward interpolation follow the
upstream formulas. The two monotone spline variants use a common
Fritsch-Carlson/Hyman-limited cubic Hermite implementation. This preserves
shape and knot values, but interior values may differ slightly from R's
specific `splinefun(method="monoH.FC")` and `splinefun(method="hyman")`
implementations.

As in R's default `approxfun`/`splinefun` setup used by the package, this port
does not extrapolate nonparametric curves.

## Nelson-Siegel fitting

R uses `optim(method="L-BFGS-B")` with analytic gradients. The Fortran port
uses a self-contained bounded Nelder-Mead method and retains the same bounds.
It minimizes the same unweighted sum of squared rate errors. Fitted parameters
can differ slightly because the optimization algorithms differ.

## Omitted infrastructure

The following are outside the computational translation:

- base-R and `ggplot2` plotting;
- S3/S4 formatting, display, data-frame, and operator machinery;
- package startup hooks and bundled datasets;
- vignettes and scripts tied to external web/data packages;
- holiday definitions from `bizdays`.
