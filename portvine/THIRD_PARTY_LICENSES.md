# Third-party components

## portvine

- Upstream project: `portvine`
- Copyright: 2021 portvine authors
- License: MIT
- Retained notice: `UPSTREAM_LICENSE.md`
- Source snapshot: `reference/portvine-master`

## rugarch Fortran translation

- Location: `vendor/rugarch`
- License: GPL-3.0-only
- Used for ARMA/GARCH types, fitting, filtering, distributions, random numbers,
  and numerical linear algebra.

## rvinecopulib Fortran translation

- Location: `vendor/rvinecopulib`
- License: GPL-3.0-only
- Used for continuous parametric pair-copulas, C-vines, D-vines, fitting,
  Rosenblatt transforms, and simulation.
- Its own third-party notices are retained under that directory.

The R package's `ppcor` dependency is not linked or copied. The small partial-
correlation calculation needed for D-vine ordering is independently implemented
in `src/portvine_ordering.f90` using standard correlation-matrix identities.
