# Testing

## FPM

```sh
fpm test
fpm run
fpm run --example regular_bond
fpm run --example day_count_comparison
```

## GNU Fortran script

```sh
./scripts/test_gfortran.sh debug
./scripts/test_gfortran.sh release
```

The debug configuration uses Fortran 2018, warnings as errors, runtime checking,
and floating-point traps. The release configuration uses optimization and the
same warning discipline.

## Test programs

- `test_daycount`: fixed independent references for all sixteen day-count
  conventions and accrued-interest calculations
- `test_schedule`: regular and long odd-coupon schedules, inferred EOM behavior,
  dates, classifications, and cash flows
- `test_pricing`: regular, odd-coupon, and zero-coupon prices; accrued interest;
  duration; convexity; and price/yield inversion
- `test_compat`: original-name exports, native date kernels, payment calculation,
  Newton solving, analytical price derivatives, duration, and convexity

## Independent references

Reference values were generated independently from the published package
formulas using Python scalar date arithmetic and NumPy/SciPy-compatible formulas.
The Fortran library does not call Python or R at runtime. See
`REFERENCE_GENERATION.md`.
