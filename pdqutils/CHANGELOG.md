# Changelog

## 0.1.0 - 2026-08-15

- Initial modern Fortran/FPM translation of PDQutils 0.1.6.
- Translated all nine exported computational APIs.
- Added standalone Edgeworth and Cornish-Fisher/AS269 implementations.
- Added generalized Gram-Charlier normal, gamma, beta, arcsine, and Wigner bases.
- Removed runtime dependencies on R, `orthopolynom`, and `moments`.
- Added scalar and vector interfaces plus Cornish-Fisher random generation.
- Corrected beta-basis integrated-polynomial shape ordering.
- Handled the arcsine degree-zero Jacobi norm analytically.
- Added robust Cornish-Fisher probability endpoint handling.
- Added strict regression, tail, vector, and RNG tests.
