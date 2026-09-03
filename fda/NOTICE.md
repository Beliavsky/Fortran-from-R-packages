# NOTICE and provenance

## Upstream

This work is derived from the computational code of the R package **fda
6.3.0**, supplied in `fda-master.zip` and identified by its `DESCRIPTION` as:

- Package: `fda`
- Version: `6.3.0`
- Date: `2025-05-21`
- Author/creator: James Ramsay
- Contributors: Giles Hooker and Spencer Graves
- Upstream URL: `http://www.functionaldata.org`
- Upstream license: `GPL (>= 2)`

The upstream `DESCRIPTION`, `NAMESPACE`, `MD5`, and `inst/NEWS.Rd` are retained
verbatim under `upstream/` for provenance.  The R source itself is not vendored
because this repository directory is a Fortran translation rather than a copy
of the R package.

## Translation license

All newly maintained Fortran source is marked `GPL-2.0-or-later`, consistent
with the upstream package declaration.  GPL version 2 is included in `COPYING`
and `LICENSE`; the "or later" permission comes from the upstream license
choice and the SPDX identifiers on translated source.

## Dependency review

Before implementation, the root of
`Beliavsky/Fortran-from-R-packages` was checked for reusable packages.
`rfortran-core`, `rfortran-linalg`, and `deSolve` were already available as
shared top-level translations/modules; no top-level `fda` or `fds` translation
was listed.

This package uses:

- `../rfortran-core` for the shared `dp = real64` kind (`r_kinds`);
- `../rfortran-linalg` for Cholesky factorization, SPD solves, inversion,
  symmetric eigensolvers, and SVD.

`rfortran-linalg` pins `fortran-lapack` at commit
`9982df3b660918761b993eaf0ff2e507547e25b2`, so no BLAS/LAPACK source is copied
here and no system `-lblas`/`-llapack` link is used.

Upstream imports `fds` and `deSolve`, but the portable numerical subset
translated here does not call those packages.  In particular, upstream
`odesolv.R` supplies its own Runge–Kutta/Cash–Karp solver, which is translated
directly.  Consequently adding unused sibling dependencies would not improve
the Fortran API and has been avoided.

## Deliberate numerical correction

The supplied upstream file `R/odesolv.R` uses `575/512` as the coefficient of
`ak2` in the state used to evaluate Cash–Karp stage 6.  The Cash–Karp tableau
uses `175/512`.  Testing the literal upstream coefficient on the elementary
problem `u'' + u = 0` caused excessive adaptive step reduction and a material
endpoint error.  The maintained Fortran therefore uses `175/512`; this is an
intentional bug fix rather than bit-for-bit reproduction of that literal.

The supplied upstream `R/monomial.R` also contains debugging `print()` calls
and its derivative recurrence multiplies by `(degree-ideriv)`, which gives an
incorrect falling-factorial coefficient (for example, the second derivative of
`x^2` is computed as zero).  It also transforms the argument without applying
the corresponding derivative scale factor.  The Fortran implementation follows
the documented mathematical derivative instead: the correct falling factorial
is multiplied by the appropriate argument-scale power.  This is another
intentional numerical correction, while the zero-derivative values retain the
upstream definition.

No other dependency source, BLAS/LAPACK/ARPACK implementation, `r.f90`,
`r_mod.f90`, or translated R-package dependency is included in this directory.
