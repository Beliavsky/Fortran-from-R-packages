# Computational coverage

This document distinguishes implemented numerical functionality from retained
but untranslated upstream functionality. The package is intentionally described
as a broad numerical-core port, not a complete replacement for every public R
interface in `actuar` 3.3-7.

## Implemented continuous distributions

| Upstream family | Fortran coverage |
|---|---|
| Pareto | density, CDF, quantile, RNG, raw moments, limited moments |
| Pareto type I | density, CDF, quantile, RNG, raw moments, limited moments |
| Burr | density, CDF, quantile, RNG, raw moments, limited moments |
| Generalized Pareto | density, CDF, quantile, RNG, raw moments, limited moments |
| Loglogistic | density, CDF, quantile, RNG, raw moments, limited moments |
| Inverse exponential | density, CDF, quantile, RNG, raw moments, limited moments |
| Inverse gamma | density, CDF, quantile, RNG, raw moments, limited moments |
| Inverse Weibull | density, CDF, quantile, RNG, raw moments, limited moments |
| Transformed gamma | density, CDF, quantile, RNG, raw moments, limited moments |
| Generalized beta | density, CDF, quantile, RNG, raw moments, limited moments |
| Gumbel | density, CDF, quantile, RNG, moments, MGF |
| Inverse Gaussian | density, CDF, quantile, RNG, moments, limited moments, MGF |

## Supplemental moments

Implemented raw moments, limited moments, or MGFs, as applicable, for:

- exponential
- normal
- beta
- gamma
- Weibull
- lognormal
- uniform
- chi-square

## Implemented discrete distributions

| Family | Coverage |
|---|---|
| Logarithmic | PMF, CDF, quantile, RNG |
| Zero-truncated Poisson | PMF, CDF, quantile, RNG |
| Zero-truncated geometric | PMF, CDF, quantile, RNG |
| Zero-truncated binomial | PMF, CDF, quantile, RNG |
| Zero-truncated negative binomial | PMF, CDF, quantile, RNG |
| Zero-modified Poisson | PMF, CDF, quantile, RNG |
| Zero-modified geometric | PMF, CDF, quantile, RNG |
| Zero-modified binomial | PMF, CDF, quantile, RNG |
| Zero-modified negative binomial | PMF, CDF, quantile, RNG |
| Zero-modified logarithmic | PMF, CDF, quantile, RNG |
| Poisson-inverse Gaussian | PMF, CDF, quantile, RNG |

## Aggregate-loss methods

Implemented:

- severity discretization from a CDF callback
  - upper discretization
  - lower discretization
  - rounding
  - unbiased first-moment adjustment
- direct discrete convolution
- exact compound summation from caller-supplied frequency probabilities
- generic Panjer `(a,b,0)` recursion
- Poisson, binomial, and negative-binomial Panjer wrappers
- compound simulation with frequency and severity callbacks
- aggregate VaR and CTE
- normal and normal-power CDF approximations
- typed aggregate distributions with mean, variance, and quantile methods

## Phase-type distributions

Implemented:

- density
- CDF
- random generation
- integer raw moments
- MGF
- native matrix exponential

The implementation accepts a transient subintensity matrix and initial
probability vector.

## Credibility

Implemented:

- Buhlmann-Straub nonparametric credibility
- Poisson-gamma conjugate credibility
- Bernoulli-beta conjugate credibility
- normal-normal conjugate credibility
- generic credibility premium combination

Not implemented:

- Hachemeister regression credibility
- general hierarchical credibility
- compound hierarchical simulation
- associated R summary and plotting infrastructure

## Grouped observations and coverage

Implemented:

- grouped mean and variance
- grouped quantiles
- linear ogive CDF
- empirical raw and limited moments
- ordinary deductible
- franchise deductible
- policy limit
- coinsurance
- inflation
- per-payment transformation

Not implemented:

- R `grouped.data` S3 class and extraction operators
- grouped histograms and graphics
- formula/data-frame adapters

## Risk theory

Implemented:

- generic adjustment coefficient from an MGF callback
- compound-Poisson adjustment coefficient
- explicit exponential-claim ruin probability
- selected Cramer-Lundberg phase-type ruin calculation

Not implemented:

- the complete `ruin` fixed-point interface
- all Sparre-Andersen interarrival models
- every discretization and interpolation option of the R interface

## Numerical support

Implemented internally:

- normal PDF/CDF/quantile
- regularized incomplete gamma functions and quantile inversion
- regularized incomplete beta and quantile inversion
- beta and gamma logarithms
- adaptive Simpson integration
- a direct Bessel-K integral
- xorshift RNG with normal, exponential, gamma, beta, Poisson, binomial,
  negative-binomial, and inverse-Gaussian generators

## Major untranslated distribution families and frameworks

The following remain in the retained upstream source but are outside this
release:

- complete Feller-Pareto distribution family
- Pareto II, III, and IV compatibility families
- inverse Burr and inverse paralogistic families
- loggamma and several transformed/inverse aliases
- minimum-distance estimation (`mde`)
- empirical-moment estimation (`emm`)
- all R vectorization, recycling, log-probability, and upper-tail conventions
- R API registration and C-callable interface
- data objects, demos, vignettes, localization, and plotting

## Source mapping

| Fortran module | Primary upstream areas |
|---|---|
| `actuar_continuous` | heavy-tailed continuous distribution R/C routines |
| `actuar_supplements` | `*Moments`, `*Supp`, limited moments, MGFs |
| `actuar_discrete` | logarithmic, zero-truncated/modified, PIG routines |
| `actuar_aggregate` | `discretize`, `exact`, `panjer`, `rcompound`, VaR/CTE |
| `actuar_phase_type` | phase-type distribution routines |
| `actuar_credibility` | `bstraub` and selected `bayes` models |
| `actuar_grouped` | grouped-data statistics and `coverage` calculations |
| `actuar_risk` | `adjCoef` and selected `ruin` calculations |
| `actuar_special`, `actuar_rng` | replacement numerical dependencies |
