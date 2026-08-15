# Changelog

## 0.1.0

- Initial modern Fortran/FPM translation of `expint` 0.2-1 computational code.
- Ported GSL-derived Chebyshev exponential-integral kernels for E1/E2.
- Ported arbitrary non-negative integer-order En evaluation.
- Added Ei via the upstream identity `Ei(x) = -E1(-x)`.
- Ported negative-argument incomplete-gamma continued fraction and recurrences.
- Added self-contained positive-argument upper incomplete-gamma evaluation.
- Added elemental APIs and explicit R-style recycling helpers.
- Added numerical tests, example, upstream provenance sources, and GPL-3-or-later license.
