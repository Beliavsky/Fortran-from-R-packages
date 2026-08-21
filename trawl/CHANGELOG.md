# Changelog

## 0.2.0

- Embedded the supplied `DEoptim-fortran` 0.1.0 translation of DEoptim 2.2-8.
- Replaced the simplified local DE/rand/1/bin optimizer with the translated
  DEoptim strategy-2 engine used by upstream `trawl`.
- Matched upstream trawl controls: `NP=10*npar`, `CR=0.5`, `F=0.8`,
  `strategy=2`, quiet tracing, and 1000 iterations by default.
- Seeded DEoptim-backed fits from the trawl RNG so `set_trawl_seed()` remains
  the package-wide reproducibility entry point.
- Added an exact wrapper-vs-direct-DEoptim fixed-seed integration test.
- Retained the complete supplied DEoptim translation under `vendor/` and added
  dependency/license provenance documentation.
- Re-ran the supplied DEoptim regression tests and all trawl tests under strict
  bounds/runtime checking and floating-point traps.

## 0.1.0

- Initial modern Fortran/FPM translation of `trawl` 0.2.2 computational code.
- Ported all non-plotting exported routines.
- Added self-contained random generation for Poisson, binomial, negative-binomial,
  and logarithmic-series laws.
- Added bounded differential evolution for the three GMM trawl fits that use
  `DEoptim` upstream.
- Added scalar all-root scanning/bisection and closed-form trawl integrals for
  intersection calculations.
- Added strict regression tests and an exponential-trawl simulation/fitting demo.
- Preserved complete GPL-3 upstream source under `upstream/`.
