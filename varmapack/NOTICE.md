# Notice

This directory is a modern Fortran translation of the computational portions
of the R package **varmapack 0.1.1**, written by Kristján Jónasson and released
under the MIT License. The upstream package metadata, selected computational
sources, R wrappers, tests, and its original third-party notice are retained
under `upstream/` for attribution and provenance.

The maintained Fortran implementation does **not** copy or compile the
upstream bundled SLICOT `sb03md-complete.F` or LAPACK-derived `rp_dpstrf.c`.
Their upstream notices are nevertheless retained verbatim in
`upstream/inst/THIRD-PARTY-NOTICES` because they are part of the provenance of
the source package from which this translation was made.

The Fortran package instead uses sibling FPM dependencies:

- `randompack` for random-number generation and multivariate-normal draws.
- `rfortran-linalg` for Cholesky factorization, symmetric eigensystems,
  spectral radii, and positive-definite linear solves.

No source from either dependency is vendored here. `rfortran-linalg` in the
target repository uses the FPM-compatible Fortran LAPACK work originating with
Federico Perini; the target repository documents its Beliavsky-maintained fork.
That dependency and its notices remain outside this package.
