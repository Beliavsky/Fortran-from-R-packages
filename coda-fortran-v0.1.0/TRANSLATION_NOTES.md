# Translation notes

## Source

Original package: `coda` 0.19-4.1.

Original package DESCRIPTION license field: `GPL (>= 2)`.

The translation follows the algorithms in the package R sources, especially:

- `R/HPDinterval.R`
- `R/autocorrdiag.R`
- `R/batchSE.R`
- `R/gelman.R`
- `R/geweke.R`
- `R/heidel.R`
- `R/mcmc.R`
- `R/mcmclist.R`
- `R/output.R`
- `R/raftery.R`
- `R/rejectionRate.R`

## Adaptation choices

R's dynamic S3 classes are represented by explicit Fortran derived types. Functions that operate on either one chain or a chain list use separate Fortran procedures where that makes the interface clearer.

Iteration labels are represented by integer `start`, `finish`, and `thin` values. This reflects the practical MCMC use in coda and avoids importing R time-series semantics.

The plotting code was deliberately not translated, as requested. Interactive menus, R printing methods, and package-specific file readers were also omitted because they are not computational algorithms.

## Compatibility caveat

Exact bit-for-bit equality with R is not promised where the R implementation delegates to `stats::ar`, `stats::glm`, `stats::qf`, or other R runtime numerical routines. The Fortran implementations use the same mathematical definitions and are tested for numerical agreement/invariants.
