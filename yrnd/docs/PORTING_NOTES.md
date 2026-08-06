# Porting notes

## Numerical calibration

The R code performs a grid search over mixture probabilities with `nlminb`, then
refines the best candidate with `constrOptim(..., method = "Nelder-Mead")`.
The Fortran port uses a deterministic multi-start bounded Nelder-Mead search over
all mixture parameters, including the component weights. This retains the same
least-squares objective, forward-mean penalty, parameter bounds, option formulas,
and American-option lower/upper-bound interpolation while avoiding thousands of
nearly identical optimizer invocations for a three-component mixture.

Mixture component labels are not economically identified. Different runs or
implementations may permute components while producing the same density and
option prices.

## Date and cash-flow representation

R `Date` objects and `lubridate` operations are replaced by `date_t` and native
Gregorian date arithmetic. The four package day-count choices are available as
`dc_act_act`, `dc_act_360`, `dc_act_365`, and `dc_30_360`.

Bond yields are solved with `tvm::xirr` from the attached modern Fortran `tvm`
translation. Coupon schedules are generated backward from maturity, supporting
annual and semiannual coupons as in the source package.

## Density construction

The R code expands a 0.001-spaced strike range iteratively until its trapezoidal
mass exceeds 0.9991. The Fortran code uses extreme component quantiles to form a
finite grid and accepts an optional `grid_step` argument. PDF values are
trapezoid-normalized; CDF values come from the analytic mixture CDF.

Moments of a lognormal mixture are evaluated analytically. Transformed STIR and
bond-yield moments are integrated numerically on their output grids.

## CTD calculations

The R implementation tests each bond as a possible CTD under a parallel shift of
all current yields and retains self-consistent candidates. The Fortran routine
implements the same rule and preserves the source fallback to the final
candidate's parallel shift when no self-consistent candidate exists.

For `proba_ctd_opt` and bond-yield transformations, coupon carry is accumulated
to the futures delivery date under the implied continuously compounded forward
rate between option expiry and futures delivery.

## Yield-spread simulation

`MASS::mvrnorm` and `stats::density` are replaced by a two-dimensional Gaussian
copula and a native Gaussian KDE. A seed can be supplied for deterministic
results.

## Intentionally omitted code

- `ggplot2` plots and formatting;
- `dplyr`, `tibble`, and `zoo` container operations;
- Bloomberg Terminal/API retrieval through `Rblpapi`;
- package documentation machinery and R-specific metadata.

The original R source remains under `original/yrnd-master` for attribution and
algorithm review.
