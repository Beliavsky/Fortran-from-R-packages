# Validation

The project is validated with GNU Fortran 14.2.0 in Fortran 2018 mode.

Checked build flags:

```text
-std=f2018 -O0 -g -Wall -Wextra -Werror -fcheck=all -fbacktrace
```

An optimized `-O2` build is also tested.

The tests cover:

- fixed independent EDGE references
- AR, AR2, CS, CS2, Roll, OHL, OHLC, CHL, and CHLO references
- combined generalized estimators
- fixed, expanding, adaptive, and endpoint windows
- multi-method dispatch
- missing-value behavior
- simulation reproducibility
- zero-observation periods
- OHLC ordering invariants
- recovery of a plausible spread from a seeded simulation

Run `scripts/validate.sh` on Unix-like systems or `scripts/validate.bat` on
Windows. FPM users can run `fpm test`.
