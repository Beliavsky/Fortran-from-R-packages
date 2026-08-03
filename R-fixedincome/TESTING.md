# Testing

Four test programs are included:

- `test_terms_compounding`: day-count conversions, term parsing, shifting,
  date offsets, compounding, and implied rates;
- `test_curves_forward`: spot factors, curve construction, forward conversion,
  round trips, helpers, and curve point insertion;
- `test_interpolation`: linear, log-linear, natural cubic, Hermite, monotone,
  and flat-forward interpolation;
- `test_parametric`: Nelson-Siegel/Svensson evaluation and synthetic fitting.

Run with FPM:

```text
fpm test
```

The GNU Fortran validation script compiles with Fortran 2018 conformance,
warning-as-error, bounds checking, uninitialized/shape checks, floating-point
traps, and backtraces:

```text
scripts/test_gfortran.sh
```

An optimized validation can be performed with:

```text
gfortran -std=f2018 -O3 -Wall -Wextra -Wpedantic -Werror ...
```
