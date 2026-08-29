# Attribution and provenance

This is a modern Fortran/FPM translation of the computational code in the R
package **statmod 1.5.2** (2026-05-17), authored by Gordon Smyth and Lizhong
Chen, with contributions from Yifang Hu, Peter Dunn, Belinda Phipson and
Yunshun Chen.

The upstream package declares `License: GPL-2 | GPL-3`. Code derived from
`statmod` therefore retains the choice **GPL-2.0-only OR GPL-3.0-only**.
Original package metadata, computational R source, native source and saved
regression-test output are retained under `upstream/`.

Important upstream numerical lineages and citations retained by this port
include:

- Dunn & Smyth (1996), randomized quantile residuals.
- Giner & Smyth (2016), inverse-Gaussian probability calculations.
- Hu & Smyth (2009), ELDA / limiting dilution analysis.
- Phipson & Smyth (2010), exact p-values for randomly sampled permutations.
- Smyth (1998/2005), Gaussian quadrature and secure-convergence nonlinear
  model fitting.
- Smyth (2002), REML for heteroscedastic regression.
- `src/04_statmod_gaussquad.f90` retains attribution to the Netlib `gaussq.f`
  / EISPACK IMTQL2 lineage described in the upstream `src/gaussq2.f`.
- `src/05_statmod_expected_deviance.f90` translates the upstream native
  Chebyshev approximation code and its coefficient tables.

`src/r_mod.F90` is the separately supplied MIT-licensed compatibility module.
The exact supplied file is retained as `upstream/r_mod-original.f90`. The
build copy differs only by whitespace/free-form continuation wrapping required
to keep source lines within the standard 132-column limit; a normalized
comparison after removal of whitespace and continuation ampersands is
identical.

`src/tweedie_dep/` is the previously translated computational core of the R
package `tweedie`, used only to provide the optional `qres.tweedie` path. It
retains its GPL-2.0-or-later license and attribution in `LICENSES/tweedie/`.
