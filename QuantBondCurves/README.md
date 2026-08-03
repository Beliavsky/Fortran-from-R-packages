# QuantBondCurves-fortran

A modern Fortran 2018/FPM translation of the computational core of the R package **QuantBondCurves 0.3.3**.

The library values fixed-income instruments and swaps, constructs coupon schedules and cash flows, converts spot and instantaneous-forward curves, and calibrates discount and cross-currency basis curves.

## Build

```sh
fpm build
fpm test
fpm run
```

GNU Fortran scripts are also provided:

```sh
./scripts/test_gfortran.sh
./scripts/test_gfortran_optimized.sh
```

On Windows:

```bat
scripts\test_gfortran.bat
scripts\test_gfortran_optimized.bat
```

## Main module

```fortran
use quant_bond_curves
```

Dates use the explicit `qbc_date` type:

```fortran
type(qbc_date) :: today
today = make_date(2026, 8, 2)
```

The core public types are:

- `qbc_date`
- `qbc_curve`
- `qbc_bond`
- `qbc_swap`
- `qbc_coupon_schedule`
- `qbc_bond_sensitivity`
- `qbc_calibration_result`

## Examples

The `example/` directory contains programs for:

- bond pricing, yield inversion, duration, convexity, and DV01;
- spot/forward curve conversion;
- bond-price curve calibration;
- interest-rate swap valuation;
- cross-currency basis-curve calibration.

## Design differences from R

- Fortran uses typed records and numeric matrices rather than data frames and mixed date/numeric matrices.
- Calendar operations are self-contained. Weekend business-day adjustment is implemented; country-specific holiday databases from `quantdates` are not bundled.
- `curve_calibration` handles the direct market-yield interpolation mode. Price-based calibration is exposed explicitly as `bootstrap_curve`.
- Calibration uses a native bounded Nelder-Mead implementation. The attached `Rsolnp` translation is GPL-2.0-only and therefore is not statically linked into this GPL-3.0-or-later FPM package.
- Plotting, vignettes, R warnings, row/column names, and R object infrastructure are not reproduced.

See `API_MAP.md`, `PORTING_NOTES.md`, and `DEPENDENCY_NOTES.md` for details.

## License

GPL-3.0-or-later. The complete upstream source snapshot is retained under `original/`.
