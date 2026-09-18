# NOTICE

This is a modern Fortran translation of computational code from the R package **kza**, version 4.2.0.

Upstream authors listed in `DESCRIPTION`:

- Brian Close
- Igor Zurbenko
- Mingzeng Sun

The upstream package declares **GPL-3** licensing. Several older compiled source files also contain their original Brian D. Close copyright notices and GNU GPL notices. Those source headers are retained verbatim in `upstream/src/`.

Important source-specific attribution retained here includes:

- `src/kz.c` upstream: KZ window functions, copyright Brian D. Close, 2015.
- `src/kza.c` upstream: adaptive KZA implementation, copyright Brian D. Close, 2005.
- `src/kzsv.c` upstream: KZ sample variance, copyright Brian D. Close, 2005.
- `R/kzft.R` upstream: copyright Brian Close, 2016, GPL version 3 notice.
- `R/rlv.R` upstream identifies Mingzeng Sun and Igor G. Zurbenko as authors of the rolling local variance implementation.

The Fortran translation is a derivative work distributed under GPL-3.0-only in this package directory. See `LICENSE` for the GNU General Public License version 3 text.

No BLAS, LAPACK, FFTW, ARPACK, `rfortran-compat`, or translated R-package dependency source is copied into this package.
