# Translation notes

## Upstream

- R package: `poibin`
- Upstream version translated: 1.6 (2024-08-23)
- Author/maintainer: Yili Hong
- Upstream license: GPL-2
- Main reference: Hong, Y. (2013), *Computational Statistics & Data Analysis*
  59, 41-51.

The upstream C source also contains the R Core Team adaptation of Singleton's
mixed-radix FFT. Upstream attribution and copyright information are retained
here and in the package metadata.

## Computational coverage

Translated:

- `dpoibin`
- `ppoibin`
  - DFT-CF exact method
  - recursive-formula exact method (`RF`)
  - refined normal approximation (`RNA`)
  - normal approximation (`NA`)
  - Poisson approximation (`PA`)
- `qpoibin`
- `rpoibin`
- integer multiplicity weights (`wts`)

R registration code and R vector/S3 infrastructure are not needed in the
Fortran package.

## Numerical implementation differences

The upstream DFT-CF path constructs the characteristic function and applies
Singleton's mixed-radix FFT. The Fortran translation constructs the same
characteristic-function samples in log-magnitude/phase form, but performs the
inverse DFT directly. This gives the same exact finite-support distribution up
to floating-point rounding and avoids carrying several hundred lines of
R-specific FFT workspace/factorization code. Its asymptotic cost is O(N^2), so
for very large support sizes the upstream FFT implementation is faster.

The `RF` method is implemented as the algebraically equivalent backward
Bernoulli recursion. Multiplicity weights are applied without first allocating
an explicitly repeated probability vector.

The Poisson approximation uses a self-contained regularized incomplete-gamma
implementation instead of R's `ppois`. Normal and refined-normal methods use
Fortran's standard `erfc` intrinsic.

Random generation follows upstream semantics: uniform draws are mapped through
the exact quantile function. The Fortran array sampler computes the exact CDF
once and reuses it for all draws.

## Invalid arguments

R signals errors for invalid probability or quantile inputs. Scalar Fortran
functions return a safe sentinel (`0` for probability functions, `-1` for an
invalid quantile) because Fortran functions do not naturally carry R-style
conditions. Low-level array constructors additionally expose integer status
codes.
