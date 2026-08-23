# Changelog

## 0.1.0

- Initial modern Fortran/FPM translation of rangen 0.0.1.
- Added portable PCG32 engine with deterministic seeding.
- Added 14 vector RNG families, matrix generators, and column-parameter generators.
- Added integer/vector/row/column sampling utilities.
- Added portable nanosecond-scaled timer.
- Replaced unavailable external Ziggurat normal backend with Box-Muller.
- Corrected upstream Cauchy, replacement-sampling, `Sample.int`, subunit-gamma, and deterministic-seeding issues.
- Added strict Fortran 2018 test suite.
