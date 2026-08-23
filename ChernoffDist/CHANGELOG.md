# Changelog

## 0.1.0

- Initial modern Fortran/FPM translation of `ChernoffDist` 0.1.0.
- Ported `dChern`, `pChern`, and `qChern`.
- Added vector wrappers and standard log/tail probability options.
- Embedded the 100 Airy zeros and derivatives needed by the upstream algorithm.
- Replaced R `integrate()` with adaptive 15-point Gauss-Kronrod quadrature.
- Replaced `uniroot()` with safeguarded Newton/bisection quantile inversion.
- Omitted unreachable coefficient entries beyond the upstream `m1=20` limit.
- Added independent numerical reference and identity tests.
