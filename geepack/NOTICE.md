# NOTICE

This directory is a modern free-form Fortran translation of the computational
parts of the R package **geepack**, version 1.3.13.

## Upstream attribution

The supplied upstream `DESCRIPTION` identifies these copyright holders and
authors:

- Søren Højsgaard — author, maintainer, copyright holder
- Ulrich Halekoh — author, copyright holder
- Jun Yan — author, copyright holder
- Claus Thorn Ekstrøm — contributor

The upstream package is licensed as **GPL (>= 3)**. This translation is
therefore distributed under **GPL-3.0-or-later**. `LICENSE` contains the GNU
General Public License version 3. The option to use later GPL versions follows
the upstream `DESCRIPTION` license declaration.

The original metadata and the computational source files used directly during
the translation are retained, unmodified, under `upstream/` for attribution and
provenance. Those files are reference material and are not compiled by FPM.

## Translation notice

The Fortran implementation is newly written code derived from the algorithms,
interfaces, formulas, and native routines in geepack. It replaces R formula and
model-frame processing, S3 methods, printing, plotting, broom/tidy helpers, and
R/C registration with typed Fortran numerical APIs.

No source code from `rfortran-core`, `rfortran-linalg`, BLAS, LAPACK, or any
translated R-package dependency is copied into this package. Those shared
components remain sibling FPM dependencies.
