# Notices and provenance

This project is a modern Fortran translation of computational routines from **wavethresh 4.7.3**, distributed by CRAN.

Upstream `DESCRIPTION` credits Guy Nason (author/maintainer) and contributors Stuart Barber, Tim Downie, Piotr Frylewicz, Arne Kovac, Todd Ogden, and Bernard Silverman. The upstream startup message states: **Copyright Guy Nason and others 1993-2022**.

The retained upstream sources contain additional source-specific attribution, including:

- `R/function.r`: multiple-wavelet routines marked **Copyright Tim Downie 1995-6** / **Copyright Tim Downie 1995-1996**.
- `src/functions.c`: code identified as originating from Arne Kovac.
- `src/functions.c`: interval-wavelet C code first written in C++ by Markus Monnerjahn and rewritten/corrected in C by Piotr Fryzlewicz while visiting the University of Bristol in 1998-9.
- Other author/source comments preserved verbatim throughout `upstream/`.

The upstream package declares `GPL (>= 2)`. This translation uses the SPDX identifier `GPL-2.0-or-later`; `LICENSE` contains GPL version 2. The supplied upstream source tree is retained under `upstream/` so its notices and attribution remain intact, except for the generated `build/partial.rdb` artifact required to be excluded from the release. Its original path and SHA-256 are recorded in `provenance/OMITTED_UPSTREAM_BUILD_ARTIFACTS.txt`.

## Shared dependencies

The translated package depends through sibling FPM paths on:

- `rfortran-core` for the shared `dp` real kind.
- `rfortran-linalg` for the shared matrix-inverse operation used by spectrum correction.
- `waveslim` for reusable 1-D and multidimensional wavelet transform steps.

No source from those dependencies is copied or vendored into this package.
