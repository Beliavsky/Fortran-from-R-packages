# NOTICE

This project is a modern Fortran translation of computational code from the R
package **dbscan**, version 1.2.6 (2026-08-24).

Upstream authors and copyright holders include Michael Hahsler and Matthew
Piekenbrock.  Sunil Arya, David Mount, and Claudia Malzer are credited upstream
as contributors.  The supplied DESCRIPTION states that the ANN library is
copyright by the University of Maryland, Sunil Arya, and David Mount, while all
other package code is copyright by Michael Hahsler and Matthew Piekenbrock.

The upstream package is licensed under GPL (>= 2).  Individual retained native
source files may contain GPL version 3 notices; those notices remain unchanged
in `upstream/src/`.  `COPYING` contains GNU GPL version 2 and this translation is
distributed under GPL-2.0-or-later, preserving all applicable upstream notices.

The upstream source bundle contains the ANN nearest-neighbor library.  ANN has
its own license and attribution; an unchanged copy of that license is retained
at `licenses/ANN-License.txt`.  ANN source is deliberately **not** copied,
vendored, compiled, or linked in the Fortran translation.  The Fortran package
uses its own exact brute-force neighbor-search implementation.

No translated dependency package, BLAS, LAPACK, ARPACK, `r.f90`, or
`r_mod.f90` is vendored in this directory.
