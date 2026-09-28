# Notices and provenance

## rARPACK upstream

This package is a modern free-form Fortran translation/compatibility
implementation based on the R package **rARPACK 0.11-0** by Yixuan Qiu, Jiali
Mei, and the ARPACK authors.

Upstream project: <https://github.com/yixuan/rARPACK>

The uploaded upstream snapshot is retained verbatim under `upstream/`. The ZIP
identifies upstream commit `a626c32e3260a81295bef29ae098268c5df9bc7c`; the
uploaded archive SHA-256 is
`03645cdebbc98eefc0aff1d33970d52aa64013caa944a608758297613a54d99e`.
The snapshot includes the CRAN `DESCRIPTION`, `NAMESPACE`, R source files,
manual pages, examples, license metadata, `AUTHORS`, `NEWS.Rd`, and historical
`COPYRIGHTS` notice. The upstream R package is licensed under
the BSD 3-Clause license.

rARPACK 0.11-0 contains no ARPACK implementation itself. Its exported R
functions `eigs`, `eigs_sym`, and `svds` directly call the corresponding
RSpectra functions. This Fortran package does not copy or translate RSpectra
source; instead it implements the documented numerical operations against the
shared `rfortran-arpack` backend already maintained in the
Fortran-from-R-packages repository.

## ARPACK backend

The FPM dependency `rfortran-arpack = { path = "../rfortran-arpack" }` supplies
explicit interfaces to the double-precision ARPACK routines `dsaupd`, `dseupd`,
`dnaupd`, and `dneupd`. No ARPACK source is copied into this package.

The historical rARPACK copyright notice records the original ARPACK copyright
as Copyright (c) 1996-2008 Rice University, developed by D.C. Sorensen,
R.B. Lehoucq, C. Yang, and K. Maschhoff, together with later arpack-ng
attributions. See `upstream/inst/COPYRIGHTS` for the preserved text.

The shared `rfortran-arpack` package carries its own BSD license and backend
provenance. Any LAPACK implementation used transitively by that shared package
is likewise not copied into this package.

## Translation implementation

New Fortran source in `src/`, tests in `test/`, and examples in `example/` are
provided under the same BSD 3-Clause terms in `LICENSE`. The implementation
uses one real kind, `dp = real64`, throughout maintained source.
