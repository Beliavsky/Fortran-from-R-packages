# NOTICE

This project is a modern Fortran translation of the numerical core of the R package `lavaan` 0.7-2 by Yves Rosseel and contributors.

Upstream package metadata declares `License: GPL (>= 2)`. This translation is therefore distributed under GPL-2.0-or-later. See `COPYING`.

Vendored dependencies:

- `pbivnorm-fortran`, GPL-2.0-or-later, based on Alan Genz's bivariate normal probability algorithm and the CRAN pbivnorm package.
- `numDeriv-fortran`, GPL-2.0-or-later, based on the CRAN numDeriv package.
- `GPArotation-fortran`, GPL-2.0-or-later, based on the CRAN GPArotation package.

Selected original lavaan R sources and package metadata are retained under `orig/` for provenance and algorithm comparison. The full original source package supplied by the user is not required at build time.
