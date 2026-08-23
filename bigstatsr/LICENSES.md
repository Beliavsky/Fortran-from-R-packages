# Licenses

## bigstatsr translation

The upstream `bigstatsr` package declares `License: GPL-3`. The translated
`bigstatsr-fortran` sources (`src/bigstatsr*.f90`) are distributed under
**GPL-3.0-only**. See `LICENSE`.

Upstream authors: Florian Prive, Michael Blum, and Hugues Aschard.
The complete supplied upstream source archive is retained under `upstream/`.

## Vendored RSpectra Fortran layer

Files under `src/vendor/rspectra*.f90` are from the modern Fortran translation
of RSpectra and are licensed under **MPL-2.0**. See `LICENSE-MPL-2.0`.
MPL-2.0 code can be distributed in this GPL-3.0 project while retaining its
file-level MPL notices.

## External numerical libraries

The default build links to system ARPACK, LAPACK, and BLAS libraries. They are
not copied into this archive. Their own licenses apply separately.
