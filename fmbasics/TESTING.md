# Testing

## FPM

```text
fpm test
fpm run
fpm run --example rates_curve_example
fpm run --example conventions_money_example
fpm run --example credit_vol_example
```

## Direct GNU Fortran validation

Unix-like systems:

```text
./scripts/test_gfortran.sh
```

Windows:

```text
scripts\test_gfortran.bat
```

The strict script uses Fortran 2018, all common warnings as errors, bounds and
runtime checking, and traps for invalid, zero-divide and overflow exceptions.

## Test programs

- `test_rates_dates`: civil dates, day counts, compounding, conversion and
  discount-factor arithmetic.
- `test_conventions`: benchmark currencies/pairs, FX dates and index dates.
- `test_curves_credit`: zero-curve loading/interpolation, survival/hazard
  conversion, CDS bootstrap and credit interpolation.
- `test_vol_money`: original volatility-surface reference values,
  multicurrency aggregation and cash-flow construction.

The rate-conversion test uses a round trip instead of the original R test's
absolute expected value because that upstream assertion specifies a tolerance
of `1e20`, making the stated constant nonbinding. The round trip tests the
actual conversion identity portably.
