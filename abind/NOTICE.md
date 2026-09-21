# NOTICE and provenance

This directory is a Fortran translation of computational code from the R package
`abind`, version 1.4-8 (package date 2024-09-08).

Upstream authors:

- Tony Plate
- Richard Heiberger

The upstream `DESCRIPTION` declares `License: MIT + file LICENSE`. The upstream
`LICENSE` identifies year 2023 and copyright holder `openaistream authors`.
That file is preserved verbatim at both `LICENSE` and `upstream/LICENSE`.
`LICENSE-MIT` contains the standard MIT license text for convenient standalone
reference; it does not replace the preserved upstream metadata.

The original R sources used for translation are retained under `upstream/R/`.
They are provenance/reference sources, not Fortran dependencies and are not
compiled by FPM.

No BLAS, LAPACK, ARPACK, `r.f90`, `r_mod.f90`, or translated dependency source is
vendored in this package. The translated implementation is dependency-free.
