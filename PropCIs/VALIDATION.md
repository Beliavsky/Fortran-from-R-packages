# Validation scope

Reference environment: Linux x86-64, gfortran and FPM from the supplied verification
toolchain; actual R with CRAN PropCIs 0.3-0. No native Windows execution is claimed.
The package has no OS-specific Fortran code or system BLAS/LAPACK link directives.

Deterministic test programs:

- `test_reference`: 255 assertions generated from unmodified PropCIs 0.3-0,
  including every computational export for which upstream returns a usable value.
- `test_bayesian_reference`: 18 independently integrated/inverted posterior
  intervals, with separate lower/upper probability checks. Includes uniform,
  Jeffreys, asymmetric, skew fractional, and all-success posterior shapes.
- `test_properties`: analytic distributions; zero/all counts; swapped samples;
  credible mass and shortest connected width; finite/NaN/infinite input behavior;
  severe-tail/shape precision-status limits; known upstream bug regressions.
- `example/intervals.f90`: eight confidence/posterior interval computations,
  including risk-ratio far below one and the formerly broken all-success Bayesian OR.

Standard runs: `fpm build`, `fpm test`, `fpm run --example`.
Checked and release runs repeat all three operations; no fast-math flags are used.
The enclosing delivery audit records exact compiler flags, stage exit codes, and
final source hashes. Build products are removed with `fpm clean --all` before
packaging. Neither R nor Python is needed for ordinary Fortran build/test/example
execution; they are needed only to regenerate the checked-in reference sources.
