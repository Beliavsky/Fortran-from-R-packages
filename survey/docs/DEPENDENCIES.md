# Dependency integration

The user supplied Fortran translations of several R dependencies. They are
preserved under `vendor/`.

## Active dependencies

### survival-fortran

Used by the survey survival wrappers. The supplied source compiles as a
standalone Fortran library; its top-level `survival.f90` facade is included in
the active vendor copy. Based on the user-provided CRAN metadata, it is marked
as the translation of `survival` 3.8-9 under LGPL (>= 2), represented here as
`LGPL-2.0-or-later`.

### minqa

Used by generic maximum pseudo-likelihood fitting. Its Powell-derived source is
compiled without warning promotion in the strict test script because the
upstream-style implementation intentionally contains exact floating-point
sentinel comparisons that trigger GNU warning diagnostics.

### numDeriv-fortran

Used for numerical Hessians in generic survey MLE fitting.

## Reference dependencies

### splines-fortran-reference

The standalone supplied splines port is retained for provenance/reference.
`survival-fortran` already contains translated modules named `splines*`; linking
both copies would create duplicate Fortran module definitions, so this second
copy is not in the active dependency graph.

### MatrixExtra-fortran

Retained as a supplied sparse-matrix reference implementation. The current
survey numerical core uses dense arrays and does not require it. Keeping it out
of the default dependency graph also avoids imposing its GPL-3-only linked
license choice on a build that uses the GPL-2-only `minqa` port.
