# Changelog

## 0.1.0

- Initial modern Fortran/FPM translation of globalOptTests 1.1.
- Translated all 50 benchmark objective functions.
- Added `go_test`, bounds, dimension, and published-optimum metadata APIs.
- Removed R `.C` dependency.
- Made Hartmann and Gulf undefined-domain source behavior safe while preserving
  intended numerical values.
- Added complete cross-language reference validation for all benchmarks.
