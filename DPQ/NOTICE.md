# Notices and attribution

This project is a modern Fortran translation of the computational portions of
R package **DPQ**, version 0.6-1, by Martin Maechler and contributors.

Upstream DPQ attribution from `DESCRIPTION` includes:

- Martin Maechler, author/maintainer and copyright holder for the main DPQ code.
- Morten Welinder, contributions involving `pgamma` code and `pdhyper`.
- Wolfgang Viechtbauer, `dtWV`.
- Ross Ihaka, `qchisq_appr` code.
- Marius Hofert, `lsum` and `lssum`.
- R Core and the R Foundation, specified native mathematical code.

The upstream `src/gamma_inc_T1006.c` is based on Algorithm 1006 by Remy
Abergel and Lionel Moisan and is separately GPL-3 licensed upstream. The
corresponding Fortran source is therefore marked GPL-3.0-or-later.

All other DPQ-derived translation sources are marked GPL-2.0-or-later, matching
upstream's per-file licensing statement. Because the distribution contains the
GPL-3 component, the combined project is distributed under GPL-3.0-or-later.

The user-supplied `r_mod.f90` remains under the MIT License. No copyright holder
was present in that supplied file, so none is invented here.

The complete upstream DESCRIPTION, LICENSE, COPYRIGHTS, R sources, native
sources, and tests are retained under `upstream/` for provenance.
