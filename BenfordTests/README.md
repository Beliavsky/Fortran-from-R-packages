# BenfordTests-fortran

Modern Fortran/FPM translation of the computational code in the R package
**BenfordTests 1.2.0** by Dieter William Joenssen and Thomas Muellerleile.

The translation preserves the package's GPL-3 license and retains the original
R and C sources under `upstream/` for provenance.

## Implemented numerical API

- first-significant-digit extraction (`significant_digit`, `significant_digits`)
- all possible leading-digit sequences (`signifd_seq`)
- Benford PMF/CDF tables (`pbenf`, `qbenf`)
- Benford random variates (`rbenf`)
- Pearson chi-square test
- Kolmogorov-Smirnov digit-frequency test
- Chebyshev-distance test
- Euclidean-distance test
- Freedman-Watson U-squared test
- Judge-Schechter mean-digit test
- JP-square correlation test
- joint-digit Hotelling/eigenvalue test, including `all`, `kaiser`, and explicit eigenvalue-index selection
- simulated null distributions for all seven simulation-based statistics
- numerical portion of `signifd.analysis`: observed frequencies, per-digit normal-approximation p-values, and confidence intervals

The graph-drawing branch of `signifd.analysis` and R S3 `htest`/printing
infrastructure are intentionally omitted.

## Build

```text
fpm build
fpm test
fpm run --example basic
```

The code is standard free-form Fortran 2018 and does not require external
numerical libraries.

## Data convention

Benford analysis is defined for nonzero finite observations. The Fortran
frequency routines ignore zero observations because an integer API has no
natural analogue of the R implementation's `NaN` significant digit for zero.
