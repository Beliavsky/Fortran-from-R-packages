# Licensing and provenance

This package translates the computational code of **optionstrat 1.4.1**,
authored and maintained upstream by **John T. Buynak <jbuynak94@gmail.com>**.
Upstream DESCRIPTION declares `License: GPL-3`. The Fortran translation and new
maintained tests/examples are distributed under **GPL-3.0-only**; the complete
GNU General Public License version 3 is included in `LICENSE`.
Upstream provides no separate copyright or NOTICE file. Authorship comments,
package metadata, R source and manual pages are retained without alteration.
No independent upstream copyright year has been inferred or added.

Official source archive:
https://cran.r-project.org/src/contrib/optionstrat_1.4.1.tar.gz

Archive SHA-256:
`5b52aee679aa76a24db07025def3a06964f3b03bc31cdf7dd1b47d8490532d6e`

CRAN publication: 2019-12-03 19:20:02 UTC. Source retrieved and translation prepared
on 2026-10-12. Metadata is in `upstream/DESCRIPTION` and `upstream/NAMESPACE`;
the exact computational source is `upstream/R/optionstrat.R`; original documentation
is in `upstream/man/`. Copies are evidence and provenance, not compiled dependencies.
The upstream MD5 inventory was checked against the downloaded archive in batch validation.

Modifications: modern free-form Fortran API and derived result types; pure elemental
operations; explicit domain checking; stable log/rate arithmetic; zero-volatility limits;
continuous-IV and standard-theta/corrected-quantile extensions; no R graphics/data frames;
no two-decimal rounding of spread analytics. README and API_COVERAGE detail source
compatibility quirks and every material interface change.

## Dependency reuse

Before implementation, the destination repository top-level package inventory and local
shared-module APIs were inspected. No existing top-level optionstrat translation or compatible
option analytics package was present. The focused `rfortran-core` API supplies `r_kinds:dp`,
`r_distributions:r_pnorm/r_qnorm/r_dnorm`, and `r_transforms:r_log1p` and is reused through the
sibling path dependency. Its separate MIT license remains with that dependency.
No rfortran-compat, translated stats package, vendor copy, BLAS, LAPACK, ARPACK,
or direct external library linking is required. Dependency source has not been copied
into this package, and shared packages have not been modified.

The original formulas implement the Black-Scholes-Merton model. The retained R source
is the direct translation source; no third-party implementation was copied into maintained
Fortran. The corrected numerical extensions are new implementations documented and tested here.
