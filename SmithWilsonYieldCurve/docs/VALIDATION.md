# Validation report

## Compiler

The package was validated with GNU Fortran 14.2.0 using Fortran 2018 mode.

Checked flags:

```text
-std=f2018 -Wall -Wextra -Werror -O0 -g -fcheck=all -fbacktrace
```

Optimized flags:

```text
-std=f2018 -Wall -Wextra -Werror -O3
```

## Test programs

1. `test_wilson_kernel`
   - zero-time identity
   - symmetry
   - three binary64 reference kernel values
2. `test_zero_coupon_fit`
   - two-instrument example from the upstream vignette
   - reference kernel weights and discount prices
   - exact calibration repricing
   - scalar and vector curve APIs
3. `test_instruments`
   - LIBOR schedule and cashflow
   - regular swap schedule and cashflows
   - off-cycle bond schedule generated backward from maturity
   - combined unique time vector, cashflow matrix, and market values
4. `test_instrument_curve`
   - end-to-end mixed LIBOR, swap, and bond fitting
   - exact repricing
   - positive discount factors and long-term UFR convergence
5. `test_errors`
   - invalid alpha
   - singular calibration system
   - unknown instrument type
   - invalid swap payment count
   - source-compatible negative-UFR warning behavior

All tests must pass in both checked and optimized builds. The final archive is
also rebuilt and retested after independent extraction.
