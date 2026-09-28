# NOTICE and provenance

## Upstream package

This project translates computational code from the R package `car`, version 3.1-5 (Date 2026-01-05; CRAN publication 2026-02-03).

Upstream title: **Companion to Applied Regression**.

Primary upstream authors are John Fox, Sanford Weisberg, and Brad Price. Upstream also credits Daniel Adler, Douglas Bates, Gabriel Baud-Bovy, Ben Bolker, Steve Ellison, David Firth, Michael Friendly, Gregor Gorjanc, Spencer Graves, Richard Heiberger, Pavel Krivitsky, Rafael Laboissiere, Martin Maechler, Georges Monette, Duncan Murdoch, Henric Nilsson, Derek Ogle, Iain Proctor, Brian Ripley, Tom Short, William Venables, Steve Walker, David Winsemius, Achim Zeileis, and R-Core.

Upstream URLs recorded by `DESCRIPTION`:

- https://github.com/bprice2652/car_repo
- https://CRAN.R-project.org/package=car
- https://z.umn.edu/carbook

The original package metadata, NAMESPACE, NEWS, R sources, and citation file used during translation are retained in `upstream/` so that numerical provenance can be audited function by function.

## License

Upstream declares `License: GPL (>= 2)`. This translation is therefore distributed under GPL-2.0-or-later. `LICENSE` contains the GNU General Public License version 2 text; the upstream declaration permits use of version 2 or, at the recipient's option, any later version.

Copyright and authorship comments present in upstream R source files are preserved verbatim in the retained files under `upstream/R/`.

## Translation provenance

The Fortran sources were produced as a source-level translation of the retained R implementation. Numerical formulas and behavior were derived from the upstream functions identified in `fpm.toml` and `docs/API_COVERAGE.md`. R plotting, interactive behavior, formula parsing, S3 object construction, and other R-specific interface machinery were not copied into the Fortran implementation.

Before adding dependencies, the target repository `https://github.com/Beliavsky/Fortran-from-R-packages` was checked for existing shared numerical modules and translated package dependencies. This package does not currently require an external Fortran dependency. In particular, it does not vendor or link BLAS, LAPACK, ARPACK, `r.f90`, `r_mod.f90`, `rfortran-compat`, or another translated R package.

The Federico Perini `fortran-lapack` implementation and the Beliavsky FPM-compatible fork are **not used by this package at present**, so no source from those projects is included.
