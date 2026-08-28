# Notices and attribution

This project is a modern Fortran translation of the computational code in
**TruncatedNormal 2.3**, authored by Zdravko Botev and Leo Belzile.  The
upstream package declares `License: GPL-3`; code in `src/` derived from it is
therefore distributed under GPL-3.0-only.  The complete GPL version 3 text is
included as `LICENSE-GPL-3.txt`.

The principal algorithms translated here are described in:

- Z. I. Botev (2017), *The normal law under linear restrictions: simulation
  and estimation via minimax tilting*, Journal of the Royal Statistical
  Society, Series B, 79(1), 125-148. DOI: 10.1111/rssb.12162.
- Z. I. Botev and P. L'Ecuyer (2015), *Efficient probability estimation and
  simulation of the truncated multivariate Student-t distribution*,
  Proceedings of the 2015 Winter Simulation Conference, 380-391.
  DOI: 10.1109/WSC.2015.7408180.

The original R and C++ computational sources, DESCRIPTION, NAMESPACE, and the
upstream tinytest file are retained in `upstream/` for provenance.

## Vendored dependency translations

The user supplied existing Fortran translations of several upstream
TruncatedNormal dependencies.  They are retained in `dependencies/` with
all of their own notices and license files:

- `nleqslv-fortran`: GPL-2.0-or-later.  This distribution uses it under GPLv3.
- `alabama-fortran`: GPL-2.0-or-later.  This distribution uses it under GPLv3.
- `spacefillr-fortran`: MIT, with its existing third-party notices.
- `qrng-fortran`: upstream choice GPL-2 or GPL-3; this distribution uses the
  GPLv3 option for qrng-derived code.
- `r_mod.F90` inside the qrng translation: separately MIT-licensed.  It is the
  previously supplied helper module with formatting-only line wrapping in the
  qrng port; see that dependency's `LICENSE-r_mod-MIT.txt` and porting notes.

Rcpp and RcppArmadillo are not dependencies of the Fortran translation.  The
small native density/Cholesky kernels that used them upstream are translated
directly to standard Fortran and BLAS/LAPACK-compatible code.
