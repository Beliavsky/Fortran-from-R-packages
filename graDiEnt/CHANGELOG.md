# Changelog

## 0.1.0

- Initial modern Fortran/FPM translation of graDiEnt 1.0.1.
- Ported all three SQG-DE adaptation strategies.
- Ported population initialization, purification, convergence checks, and
  trace bookkeeping.
- Added explicit typed Fortran objective callbacks.
- Added deterministic Fortran RNG seeding.
- Added strict tests and examples.
- Added defensive fixes for upstream parent-count, thinning, convergence, and
  zero-jitter edge cases.
