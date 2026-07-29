# Porting notes

## Preserved algorithms

The translation follows the executable R source in `R/`:

- moment-based starting values used by `gev`, `gumbel`, `gpd`, and `pot`;
- GEV, Gumbel, GPD, and point-process negative log likelihoods;
- GPD probability-weighted-moment estimator and expected covariance formula;
- logistic bivariate POT likelihood, locally or globally fitted margins;
- source definitions of Hill, extremal-index, record, threshold, mean-excess,
  shape-stability, quantile-stability, return-level, and risk-measure routines;
- source block construction, including a final partial block;
- source threshold convention: the next distinct observation below the
  requested upper order statistic.

## Numerical substitutions

R's `optim` and `solve` are replaced by a self-contained Nelder-Mead optimizer,
central finite-difference Hessian, and pivoted Gauss-Jordan inversion. This can
produce slightly different iteration counts, convergence flags, covariance
estimates, and last digits while targeting the same likelihood.

The Fortran code handles `xi=0` by the exact Gumbel or exponential limit. The R
source divides by `xi` directly in several places, so near-zero fits can be
unstable there.

## Corrections and clarifications

1. The Wald covariance term in R `gpd.q` multiplies the parameter covariance by
   the beta variance. The translated delta method uses the mathematically
   consistent cross term
   `2 * cov(xi,beta) * g * beta * g_prime`.
2. GEV/GPD distribution functions explicitly enforce support and scale
   constraints and implement finite endpoint behavior.
3. Point-process fitting rejects nonpositive observation spans rather than
   permitting divisions by zero.
4. Bivariate fitting validates all log-likelihood factors before taking logs.
5. Expected shortfall is reported as infinite for `xi >= 1`, where the GPD
   first moment does not exist.

## R-only features omitted

- interactive plot menus and graphics;
- S3 classes and plot dispatch;
- POSIX/date parsing and `ts` attributes;
- console warnings, labels, and formatted printing;
- automatic spline drawing for profile likelihoods.

Profile grids and diagnostic coordinates are returned numerically, so callers
can plot them with any Fortran or external graphics system. Calendar blocks can
be prepared as explicit integer group identifiers and passed to
`block_maxima_groups`.
