# Notices and attribution

This directory is a modern free-form Fortran translation of computational
algorithms from the R package **strucchange**, version 1.6-0, "Testing,
Monitoring, and Dating Structural Changes".

## Upstream authors and contributors

The upstream `DESCRIPTION` identifies these authors:

- Achim Zeileis
- Friedrich Leisch
- Kurt Hornik
- Christian Kleiber

and these contributors:

- Bruce E. Hansen
- Edgar C. Merkle
- Nikolaus Umlauf

The original package metadata is retained in `upstream/DESCRIPTION`, and its
citation information is retained in `upstream/CITATION`.

## License

The upstream package declares `GPL-2 | GPL-3`. This translation preserves that
license choice. Full copies of GNU GPL version 2 and GNU GPL version 3 are in
`LICENSE-GPL-2` and `LICENSE-GPL-3`.

## Computational provenance

The Fortran code is based on computational material in the following upstream
files:

- `src/strucchange_functions.c`: recursive residual calculation.
- `R/recresid.R`: recursive residual reference implementation and indexing.
- `R/Fstats.R`: sequence of Chow/F statistics and test aggregations.
- `R/breakpoints.R`: RSS triangle, dynamic programming, BIC, segmented fits,
  breakpoint confidence intervals, and `pargmaxV`.
- `R/matrix.R`: matrix square-root and cross-product helpers.
- `R/efp.R`: OLS/recursive CUSUM, OLS/recursive MOSUM, recursive/moving
  estimates, score fluctuation processes, process functionals, and asymptotic
  fluctuation-process p-values.
- `R/gefp.R`: generalized score/estimating-function fluctuation process.
- `R/pvalue.Fstats.R`: Hansen-style response-surface coefficients used for
  supF/aveF/expF asymptotic p-values.
- `R/critvals.R`: Brownian-motion/bridge critical-value tables.
- `R/critvals-monitoring.R`: monitoring critical-value tables.
- `R/monitoring.R`: monitoring boundaries and critical-value formulas.
- `R/zzz.R`: reusable `supLM`, maximum-MOSUM, and categorical L2 functionals.

Large numerical tables in `src/strucchange_tables.f90` are transcriptions of
the corresponding upstream R data objects. They are data/provenance-bearing
content, not independently generated approximations.

No source code from `rfortran-core`, `rfortran-linalg`, BLAS, LAPACK, ARPACK,
or another translated R package is copied into this package. Those facilities
are referenced only through FPM dependencies where needed.
