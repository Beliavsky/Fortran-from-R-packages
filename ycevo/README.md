# ycevo-fortran

A modern Fortran 2018 translation of the computational code in the R package
`ycevo` 0.2.1.9000, **Nonparametric Estimation of the Yield Curve Evolution**.
The original package implements the estimator described by Koo, La Vecchia,
and Linton (2021).

The port retains the GPL-3 license and includes the original computational R
and Rcpp source files under `original/` for attribution and comparison.
Plotting, `ggplot2`, tidy-data manipulation, S3 method dispatch, and parallel
progress reporting are not translated.

## Implemented numerical functionality

- Epanechnikov density, quantile, random generation, and kernel weights.
- Date, maturity, and optional covariate kernel weighting.
- Nelson-Siegel yield curves with cubic time evolution.
- Flattened bond-panel representation replacing lists of sparse R matrices.
- `dbar` numerator and denominator calculation.
- Coupon cash-flow cross-product (`H-hat`) calculation.
- Interpolation of `H-hat` and solution of the discount-function system.
- Single-curve and multi-date yield-curve estimation.
- Optional one-dimensional covariate weighting, such as a short rate.
- Automatic default maturity grids and bandwidths.
- Counting and removal of maturity windows with too few maturing bonds.
- Local-quadratic LOESS-like smoothing and clamped linear/bilinear prediction.
- Synthetic coupon-bond panel generation.
- Bond-panel CSV input and yield-curve CSV output.

## Build with FPM

```text
fpm test
fpm run --example ycevo_example
```

The example is declared as an executable, so with current FPM releases use:

```text
fpm run ycevo_example
```

## Build without FPM

On Unix-like systems:

```text
./scripts/test_gfortran.sh
```

For an optimized build:

```text
./scripts/build_gfortran.sh
./build/bin/test_ycevo
./build/bin/ycevo_example
```

## Core API example

```fortran
use ycevo, only : dp, ycevo_success
use ycevo, only : bond_panel_t, yield_curve_t
use ycevo, only : simulate_bond_panel, estimate_yield

type(bond_panel_t) :: bonds
type(yield_curve_t) :: curve
real(dp) :: tau(4), ht(4)
integer :: status
character(len=256) :: message

call simulate_bond_panel(bonds, nday=40, n_bonds=60, seed=12345)
tau = [0.5_dp, 1.0_dp, 2.0_dp, 5.0_dp]
ht  = [0.4_dp, 0.5_dp, 0.7_dp, 1.0_dp]
call estimate_yield(bonds, 0.5_dp, 0.35_dp, tau, ht, &
                    curve, status, message)
if (status /= ycevo_success) error stop trim(message)
```

## Bond-panel layout

`bond_panel_t` contains one row per future cash flow:

- `day`: integer quotation-day index in `1:nday`.
- `id`: integer bond identifier.
- `price`: bond price, repeated for each cash-flow row of a bond/day.
- `tupq`: days from quotation to payment.
- `cashflow`: coupon or principal-plus-coupon payment.

The CSV reader expects these columns in this order:

```text
day,id,price,tupq,cashflow
```

## References

- Koo, B., La Vecchia, D., and Linton, O. B. (2021). Estimation of a
  nonparametric model for bond prices from cross-section and time series
  information. *Journal of Econometrics*, 220(2), 562-588.
- Koo, B., and Yang, Y. F. (2024). `ycevo`: Non-Parametric Estimation of the
  Yield Curve Evolution, R package version 0.2.1.

See `docs/API_MAP.md` and `docs/PORTING_NOTES.md` for translation details.
