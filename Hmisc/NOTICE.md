# NOTICE and provenance

This directory is an independent modern Fortran translation of selected computational code from **Hmisc 5.3-0** (upstream package date 2026-09-05), authored principally by Frank E. Harrell Jr with contributors including Cole Beck and Charles Dupont.

Upstream source: <https://github.com/harrelfe/Hmisc> and CRAN package Hmisc.

The upstream package is licensed GPL (>= 2).  The translated code is distributed under GPL-2.0-or-later; see `LICENSES/GPL-2.0-or-later.txt`.  Algorithmic provenance is retained in `fpm.toml` and `README.md` by mapping translated R functions to their upstream source files.  Several routines (`rcorr`, `hoeffd`, `rcorr.cens`, `rcorrp.cens`, `cutGn`, `whichClosest`, and `whichClosePW`) correspond to public R wrappers around compiled kernels shipped by Hmisc; this translation re-expresses those numerical operations in a module-oriented Fortran API.

No BLAS, LAPACK, ARPACK, r.f90, r_mod.f90, translated dependency package, or other vendored numerical dependency is included.  The current translated subset does not require one.
