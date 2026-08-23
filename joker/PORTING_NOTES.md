# Porting notes

## Upstream source

Translated from `joker` 0.14.2. The original DESCRIPTION, NAMESPACE and R
sources are retained under `orig/` for provenance.

## Numerical substitutions

R's `stats` distribution primitives are replaced by self-contained Fortran
implementations using regularized incomplete beta/gamma functions, bisection
quantiles, and native random generators.

The upstream beta MLE reduces optimization to the total concentration and then
uses inverse-digamma updates. The Fortran implementation solves the equivalent
two-parameter score equations with Newton iteration.

The upstream gamma and Weibull MLEs are implemented through their one-dimensional
score equations. Dirichlet MLE uses the standard Newton update for the complete
Dirichlet score. Multivariate-gamma MLE uses the independent gamma increments
implied by the cumulative representation and solves the shared-scale score
system.

## Numerical parity added in v0.2.0

The complete collection of upstream estimator asymptotic covariance methods is
now translated: every `avar_mle`, `avar_me`, and `avar_same` method present in
the R package has a Fortran counterpart. A new `joker_asymptotics` module also
exposes the corresponding Fisher-information helpers.

Three clear source-level information-matrix errors are corrected rather than
reproduced:

- `Multinom::finf` uses a minus sign for the constrained-simplex term. The
  Fisher information requires a plus sign, yielding
  `Cov(p_hat) = (diag(p) - p p^T) / size` for the first `k-1` probabilities.
- `Laplace::finf` omits one power of `sigma`. The correct information for both
  location and scale is `1 / sigma^2`.
- `Lnorm::finf` divides by `sdlog` rather than `sdlog^2`. Since log-normal
  likelihood in `(meanlog, sdlog)` is the normal likelihood of `log(x)`, the
  information is `diag(1, 2) / sdlog^2`.

The upstream tests use unit scale for the latter two cases, which masks those
errors. The v0.2 tests explicitly use non-unit scales.

The R implementation applies `Matrix::nearPD` to several covariance matrices.
The Fortran port returns the analytical covariance formula directly and avoids
a Matrix dependency; for regular interior parameter values these matrices are
positive semidefinite up to floating-point roundoff.

## Intentionally omitted after v0.2.0

- S4 class dispatch and reflection.
- plotting and ggplot2/ggh4x presentation.
- `SmallMetrics` / `LargeMetrics` data-frame simulation/reporting wrappers.
- R input-validation helpers and data-frame presentation wrappers.
- presentation-only `moments()` aggregators.

Fisher-F, Student-t, uniform, and Weibull do not define `avar_*` methods in
upstream joker 0.14.2, so their absence is not a covariance-parity gap.
