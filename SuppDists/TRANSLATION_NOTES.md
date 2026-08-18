# Translation notes

## Upstream

This project translates the computational code of CRAN package SuppDists
1.1-9.9 (2025-03-23), authored by Bob Wheeler with contributions and current
maintenance by Thorsten Pohlert. The upstream package declares `GPL (>= 2)`;
this translation is distributed under GPL-2.0-or-later. The upstream
DESCRIPTION is retained in `licenses/UPSTREAM_DESCRIPTION`.

The original implementation is primarily `src/dists.cpp`, `src/dists.h`,
`src/datatabs.h`, and the R wrappers in `R/Dists.R`.

## What was translated

The exported computational families are represented by Fortran procedures:

- `dinvgauss`, `pinvgauss`, `qinvgauss`, `rinvgauss`, `sinvgauss`
- `dkruskalwallis`, `pkruskalwallis`, `qkruskalwallis`, `rkruskalwallis`,
  `skruskalwallis`
- `dnormscore`, `pnormscore`, `qnormscore`, `rnormscore`, `snormscore`,
  `norm_order`
- `dkendall`, `pkendall`, `qkendall`, `rkendall`, `skendall`
- `dfriedman`, `pfriedman`, `qfriedman`, `rfriedman`, `sfriedman`
- `dspearman`, `pspearman`, `qspearman`, `rspearman`, `sspearman`
- `dmaxfratio`, `pmaxfratio`, `qmaxfratio`, `rmaxfratio`, `smaxfratio`
- `dpearson`, `ppearson`, `qpearson`, `rpearson`, `spearson`
- Johnson SN/SL/SU/SB density/CDF/quantile/RNG/statistics plus
  `johnson_fit_quantiles` and `johnson_fit_moments`
- `dghyper`, `pghyper`, `qghyper`, `rghyper`, `sghyper`, `hyper_type`, and
  `hyper_type_name`
- `moments` and `sample_moments`

The exact Friedman/Spearman tables in upstream `datatabs.h` were converted to
Fortran and are used for the same tabulated ranges: Friedman r=3 through the
upstream n limits, r=4 and r=5 through their upstream limits, and the n=2
tables through r=11 used for exact Spearman probabilities. Outside those
ranges the same beta-family approximation used upstream is applied.

Kendall uses the upstream inversion-count dynamic recursion exactly for
N <= 12 and the upstream Edgeworth expansion for larger N.

Kruskal-Wallis and normal-score distributions use the Wallace and Lu-Smith
variance formulas and the same beta approximation as upstream. Expected
normal order statistics translate AS 177.3 and its finite-sample corrections.

Pearson correlation uses the same Johnson-Kotz hypergeometric-series density.
Its CDF is numerically integrated and the quantile is inverted by bisection.
The maximum-F CDF and density use Hartley's defining chi-square integrals.

## Deliberate implementation differences

The R `.C` vector-recycling layer, S3/list formatting, `lower.tail`, `log.p`,
and `log` wrapper arguments are not reproduced as language infrastructure.
The Fortran API is scalar-first; callers can loop or use array wrappers of
their own. Quantile procedures are lower-tail quantiles.

The generalized hypergeometric implementation identifies the same Kemp-Kemp
types but normalizes the probability recurrence directly. This avoids several
separate gamma-function branches in the C++ implementation. Infinite-support
families are summed until the next term is negligible (with a large safety
cap).

Johnson quantile fitting follows Wheeler's quantile construction. Normal and
lognormal moment fits use the upstream closed forms; SU/SB moment fitting uses
numerical matching of standardized skewness and kurtosis rather than a
line-for-line translation of AS 99. This preserves the target four moments
but may produce slightly different fitted parameters when several Johnson
curves fit almost equally well.

Numerical integration uses adaptive Simpson quadrature rather than the
upstream integration helper. Quantile inversion generally uses robust
bisection rather than Newton iteration.

Random generation uses Fortran's `random_number` plus Box-Muller / standard
transformation algorithms. The old MWC1019 and Ziggurat code is commented out
and unregistered in the supplied upstream source, so it is not exposed here.

## Tests

The test suite contains independent numerical checks for inverse Gaussian,
classic and generalized hypergeometric probabilities, exact Kendall,
upstream Friedman/Spearman tables, Kruskal-Wallis, Johnson transformations,
Pearson correlation, maximum F ratio, and normal-order-score symmetry.
