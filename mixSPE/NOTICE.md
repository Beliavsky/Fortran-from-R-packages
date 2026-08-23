# Licensing and provenance

This is a modern Fortran translation of the computational code in the R package
`mixSPE` 0.9.3 (2025-04-08), by Ryan P. Browne, Utkarsh J. Dang,
Michael P. B. Gallaugher, and Paul D. McNicholas.

The upstream package declares `License: GPL (>= 2)`.

The user supplied `mvtnorm-fortran` translation is used for the multivariate
normal proposal distribution required by `rspe`. Its source modules used here
retain their SPDX headers, and its GPL-2.0-only license is copied under
`vendor/mvtnorm/LICENSE`. Consequently this combined source distribution is
provided under GPL-2.0-only.

No plotting or R S3/data-frame presentation code is translated.
