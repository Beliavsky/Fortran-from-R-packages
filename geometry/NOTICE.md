# NOTICE and provenance

This package is a modern free-form Fortran translation of computational code from the R package `geometry` version 0.5.2.

Upstream package metadata names Jean-Romain Roussel, C. B. Barber, Kai Habel, Raoul Grasman, Robert B. Gramacy, Pavlo Mozharovskyi, and David C. Sterratt as copyright holders/authors/contributors. The upstream R package is distributed under GPL (>= 3). This translation is distributed under GPL-3.0-or-later; see `LICENSE`.

The original R sources used for translation are retained under `upstream/R/`, with their original copyright and attribution comments. In particular, coordinate-conversion routines were adapted upstream from GNU Octave code by Kai Habel; mesh-distance and DistMesh routines credit Raoul Grasman and the MATLAB work of Per-Olof Persson; `polyarea` credits David M. Doolin and the Octave source; `tsearch`/barycentric code credits David Bateman and David Sterratt; and `surf.tri` credits Per-Olof Persson and Raoul Grasman. See the retained upstream source files for the complete notices.

The upstream R distribution embeds the Qhull C library and documents its separate notices in `upstream/LICENSE-NOTES`. This Fortran translation does **not** copy, vendor, translate, or link that Qhull source. The convex-hull, Delaunay, halfspace, and search routines in `src/` are independent Fortran implementations written for this translation. Qhull remains important upstream provenance because several R interfaces and their documented semantics originated as Qhull wrappers.

No source from BLAS, LAPACK, ARPACK, `magic`, `lpSolve`, `linprog`, Rcpp, or RcppProgress is included in this package.
