# sn-fortran

Modern Fortran 2018 translation of the computational core of the R package
`sn` 2.1.3 by Adelchi Azzalini. The library covers the skew-normal (SN),
extended skew-normal (ESN), skew-t (ST), skew-Cauchy (SC), multivariate SN/ST/SC,
and unified skew-normal (SUN) families.

The project is an FPM package and is also buildable with the included Makefile.
It is licensed under GPL-2.0-only OR GPL-3.0-only, matching the upstream choice.

## Main capabilities

- Univariate SN/ESN, ST, and SC density, CDF, quantile, and random generation
- Owen's T function, zeta derivatives, cumulants, modes, and quantile summaries
- Direct, centered, and operational parameter conversions
- Multivariate SN/ST/SC density, CDF, random generation, moments, marginals,
  affine transforms, and SN conditioning
- SUN density, CDF, random generation, moments, marginals, affine transforms,
  convolution, joins, equality conditioning, inequality conditioning, and
  SN/CSN conversion
- Regression likelihood fitting for normal, SN, ESN, ST, and SC errors
- Q-penalized or unpenalized fits, weighted observations, numerical covariance,
  fitted values, residuals, and prediction
- Grouped-data likelihood fitting
- Composite multivariate SN/ST fitting
- Symmetric-modulated distributions
- Products of correlated normal or Student-t variables
- `vech`, inverse-`vech`, duplication matrices, traces, and block diagonals
- Self-contained probability, random-number, optimization, and linear-algebra code

## Build

With FPM:

```text
fpm test
fpm run
```

With GNU Make:

```text
make check
make optimized
make demo
```

The checked build uses Fortran 2018 conformance, all common warnings as errors,
bounds/runtime checking, and backtraces. The optimized build uses `-O3`.

## Minimal example

```fortran
program example
  use sn, only : dp, dsn, psn, qsn
  implicit none
  real(dp) :: p

  p = psn(0.5_dp, alpha=2.0_dp)
  print *, dsn(0.5_dp, alpha=2.0_dp)
  print *, qsn(p, alpha=2.0_dp)
end program example
```

See `app/demo_sn.f90` for regression fitting and random generation.

## Numerical scope

The closed-form and low-dimensional formulas are direct translations. The
upstream package delegates higher-dimensional normal and Student-t probabilities
to `mnormt`; this port instead uses deterministic antithetic Halton integration.
SUN moments use deterministic truncated-normal simulation. Results are therefore
reproducible but are not expected to be bit-for-bit identical to R or `mnormt`.

R formulas, S3/S4 classes, plotting, model frames, data frames, and interactive
profiling are intentionally omitted. Typed Fortran parameter/result structures
replace R objects. See `docs/PORTING_NOTES.md` and `docs/API_MAP.md`.

## Provenance

The unmodified source archive and extracted upstream package are retained under
`upstream/`. Copyright remains with the original authors and contributors.
