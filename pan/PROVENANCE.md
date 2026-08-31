# Provenance

## Source package

* Package: `pan`
* Version: 2.0
* Date: 2026-06-30
* Repository field: CRAN
* License field: GPL-3

The user supplied archive was `pan-master.zip`.

Primary computational sources inspected:

* `src/pan.f` - legacy fixed-form Fortran implementation of the Gibbs samplers,
  dense matrix helpers, random-number generators, and ECME routines.
* `R/pan.R` - R wrappers for `pan`, `pan.bd`, and `ecme`, including missingness
  pattern construction and prior packing.
* `man/pan.Rd`, `man/pan.bd.Rd`, and `man/ecme.Rd` - documented statistical
  model, prior interpretation, outputs, and assumptions.

SHA-256 of the supplied primary sources:

* `src/pan.f`: `3efbf063227605275a95378627beb787bb556fa31ed367d1df819b1ecd6d45f3`
* `R/pan.R`: `a1eeefe8c5d75239ee28f225127cb0d6eb77c2f67f24a66ae55ab86138271932`

`UPSTREAM_DESCRIPTION` and `UPSTREAM_CITATION` are retained from the supplied
archive.

## Repository dependency check

Before implementation, the current root of
`Beliavsky/Fortran-from-R-packages` was checked. The repository documents
shared `rfortran-core` and `rfortran-linalg` packages, including a pinned
pure-Fortran LAPACK backend for shared linear algebra.

This translation does not require either package: the upstream `pan.f`
computational source itself includes the small dense Cholesky, inversion,
Wishart, and RNG kernels needed by its algorithms. Translating those
package-specific kernels directly avoids introducing a dependency or copying
shared-package source.

No BLAS, LAPACK, ARPACK, `r.f90`, `r_mod.f90`, or translated dependency source
is included.

## Translation approach

The legacy fixed-form entry points and workspace arrays were not copied into
the maintained source tree. Their statistical operations were reorganized into
free-form modules with allocatable result/state types and explicit interfaces.

The translation preserves the full and block-diagonal Gibbs models and the
Gaussian mixed-model ML target. See `API_COVERAGE.md` for algorithm-level
differences that do not alter the model being fitted.
