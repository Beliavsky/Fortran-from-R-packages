# Porting notes

## Upstream

- Package: AdequacyModel
- Version: 2.0.0
- License: GPL (>= 2), represented here as GPL-2.0-or-later
- Authors: Pedro Rafael Diniz Marinho, Marcelo Bourguignon, Cicero Rafael Barros Dias
- CRAN publication date in the supplied source: 2016-05-20

## API translation

R function objects are represented by Fortran procedure callbacks. R lists
returned by optimizers and `goodness.fit()` are represented by derived types.
R formula/S3 machinery is not involved in this package.

## Deliberate PSO correction

The upstream `pso.R` does not perform a standard per-particle personal-best
update. In particular, the conditional around `f_xi <= f_pi` is ineffective,
and the following block updates only the particle at a global minimum index.
It also compares values to `NaN` directly, which is not a reliable NaN test.

The Fortran routine implements the intended standard PSO semantics:

- every particle updates its own personal best when its objective improves;
- the global best is updated from all personal bests;
- IEEE non-finite objective values are handled explicitly;
- out-of-bound coordinates are resampled within their own parameter bounds;
- random cognitive/social coefficients are generated per particle and
  coordinate;
- a finite `max_iter` safeguard is available.

The upstream constants `omega=0.5`, `phi_p=0.5`, and `phi_g=0.5` and its
variance-of-recent-best stopping idea are retained.

## `goodness.fit()`

The transformed-normal Cramer-von Mises and Anderson-Darling calculations are
preserved. The Fortran code clamps CDF values away from exactly zero and one
before applying the inverse normal and logarithms. This prevents infinities at
both tails; upstream only replaces positive infinity in one intermediate
vector.

The one-sample KS statistic is computed directly against the supplied CDF.
For moderate sample sizes and matrix dimensions the p-value uses the exact
Marsaglia-Durbin matrix calculation; larger cases use the usual corrected
Kolmogorov asymptotic approximation.

PDF normalization checks support finite, semi-infinite, and doubly-infinite
domains using adaptive Simpson integration after a change of variables.

The `CAIC`/`CAIC ` naming inconsistency in the R return list is normalized to
`aicc`, since the formula in the package is the small-sample corrected AIC,
not the commonly named consistent AIC.

## Local optimization methods

The package delegates BFGS, Nelder-Mead, CG, and SANN to base R `optim()`.
This port supplies native Fortran implementations so the library is
self-contained. They optimize the same objective but are not intended to
reproduce R `optim()` iteration-for-iteration.

## Descriptive mode

For integer-valued samples, all tied exact modes are returned and the behavior
matches the upstream intent. For non-integer data, upstream obtains the mode
from the midpoint(s) of the tallest default R `hist()` bin. The Fortran port
uses a Sturges equal-width histogram without R's graphics-oriented `pretty()`
break adjustment. Therefore a continuous-data mode can differ slightly.

## Reentrancy

`goodness_fit()` temporarily stores the supplied PDF procedure pointer so it
can be passed through the generic optimizer callback interface. Calls to
`goodness_fit()` should not be nested or invoked concurrently from multiple
threads. `goodness_from_mle()` has no such limitation.
