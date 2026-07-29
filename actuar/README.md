# actuar-fortran

A dependency-free modern Fortran/FPM numerical-core port of selected reusable
algorithms from the R package `actuar` 3.3-7.

The port focuses on scalar probability laws, aggregate-loss calculations,
phase-type distributions, credibility, grouped-data calculations, insurance
coverage transformations, and selected risk-theory routines. It is not a
mechanical translation of every R, C, documentation, plotting, data, or fitting
facility in the upstream package.

## Requirements

- A Fortran 2018 compiler
- FPM for the normal build workflow
- No BLAS, LAPACK, statistics, optimization, or special-function library

The validation scripts were exercised with GNU Fortran 14.2.0.

## Build and test

```text
fpm build
fpm test
fpm run actuar_demo
fpm run --example loss_distributions
fpm run --example aggregate_and_credibility
```

Direct validation scripts are also provided:

```text
tools/validate.sh
```

or on Windows:

```text
tools\validate.bat
```

## Main modules

The simplest import is:

```fortran
use actuar
```

The umbrella module re-exports the public interfaces of these modules:

- `actuar_continuous`: heavy-tailed and transformed continuous laws
- `actuar_supplements`: moments, limited moments, and MGFs for common laws
- `actuar_discrete`: logarithmic, zero-truncated, zero-modified, and PIG laws
- `actuar_aggregate`: discretization, convolution, Panjer recursion, and risk
- `actuar_phase_type`: phase-type probabilities, moments, and simulation
- `actuar_credibility`: Buhlmann-Straub and conjugate Bayesian credibility
- `actuar_grouped`: grouped observations and insurance coverage transforms
- `actuar_risk`: adjustment coefficients and selected ruin probabilities
- `actuar_special`, `actuar_rng`: native numerical support

All floating-point calculations use:

```fortran
integer, parameter :: dp = kind(1.0d0)
```

## Continuous loss distributions

The port includes density, distribution, quantile, random generation, moments,
and limited moments where mathematically available for:

- Pareto/Lomax
- Pareto type I
- Burr
- generalized Pareto/beta-prime
- loglogistic
- inverse exponential
- inverse gamma
- inverse Weibull
- transformed gamma
- generalized beta
- Gumbel
- inverse Gaussian

Examples:

```fortran
use actuar, only : dp, dpareto, ppareto, qpareto, levpareto

real(dp) :: density, probability, quantile, limited_mean

density = dpareto(2.0_dp, shape=3.0_dp, scale=4.0_dp)
probability = ppareto(2.0_dp, shape=3.0_dp, scale=4.0_dp)
quantile = qpareto(0.95_dp, shape=3.0_dp, scale=4.0_dp)
limited_mean = levpareto(10.0_dp, order=1.0_dp, shape=3.0_dp, scale=4.0_dp)
```

Supplementary actuarial moments are available for exponential, normal, beta,
gamma, Weibull, lognormal, uniform, and chi-square distributions.

## Frequency distributions

The following scalar PMF/CDF/quantile/random interfaces are included:

- logarithmic
- zero-truncated Poisson, geometric, binomial, and negative binomial
- zero-modified Poisson, geometric, binomial, negative binomial, and logarithmic
- Poisson-inverse Gaussian

Example:

```fortran
use actuar, only : dp, dztpois, pztpois, qztpois

print *, dztpois(3, 2.5_dp)
print *, pztpois(3, 2.5_dp)
print *, qztpois(0.95_dp, 2.5_dp)
```

## Aggregate-loss calculations

Severity distributions can be discretized from a caller-supplied CDF. The
resulting probability vector can be combined using direct convolution, exact
compound summation, or Panjer recursion for Poisson, binomial, and negative
binomial frequencies.

```fortran
use actuar

real(dp), allocatable :: severity(:)
type(aggregate_distribution) :: aggregate

severity = [0.20_dp, 0.50_dp, 0.30_dp]
aggregate = panjer_poisson(severity, lambda=2.0_dp, step=100.0_dp, max_n=40)

if (.not. aggregate%ok) error stop trim(aggregate%message)

print *, aggregate%mean()
print *, aggregate%variance()
print *, aggregate%quantile(0.99_dp)
print *, aggregate_cte(aggregate, 0.99_dp)
```

Also included are compound simulation callbacks, aggregate VaR/CTE, and normal
and normal-power approximations.

## Phase-type distributions

The module supports transient generator matrices and initial probabilities:

```fortran
use actuar

real(dp) :: alpha(2), subintensity(2,2)

alpha = [1.0_dp, 0.0_dp]
subintensity = reshape([-2.0_dp, 0.0_dp, 2.0_dp, -3.0_dp], [2,2])

print *, dphtype(1.0_dp, alpha, subintensity)
print *, pphtype(1.0_dp, alpha, subintensity)
print *, mphtype(1, alpha, subintensity)
```

A native scaling-and-squaring matrix exponential is used; no external linear
algebra library is required.

## Credibility

Buhlmann-Straub credibility returns typed means, credibility weights, estimates,
and process/structural variance estimates. Conjugate Bayesian credibility is
provided for Poisson-gamma, Bernoulli-beta, and normal-normal models.

```fortran
use actuar

type(credibility_result) :: fit
real(dp) :: claims(3,4), weights(3,4)

fit = buhlmann_straub(claims, weights)
if (.not. fit%ok) error stop trim(fit%message)
print *, fit%estimates
```

## Grouped data and policy coverage

Grouped means, variances, quantiles, ogive values, empirical moments, limited
moments, and common payment transformations are included. `apply_coverage`
supports ordinary or franchise deductibles, policy limits, coinsurance, and
inflation.

## Risk theory

The risk module includes:

- generic adjustment-coefficient root solving from a severity MGF callback
- a compound-Poisson adjustment-coefficient wrapper
- the explicit exponential-claim Cramer-Lundberg ruin probability
- a selected phase-type Cramer-Lundberg ruin calculation

The phase-type ruin routine should be treated as a focused implementation, not
as a replacement for every Sparre-Andersen and renewal-equation method exposed
by upstream `ruin`.

## Numerical conventions

- APIs are scalar unless an array is explicitly part of the model.
- Probabilities are lower-tail, ordinary-scale values.
- Invalid parameters generally return IEEE NaN or a typed result with
  `ok=.false.` rather than raising an R condition.
- Random generators support deterministic seeding through `seed_rng` or the
  higher-level seed arguments.
- Boundary comparisons are explicit, which is why direct validation suppresses
  GNU Fortran's warning for intentional real equality tests.

## Deliberately excluded

The current release does not claim complete coverage of the upstream package.
Notable exclusions include:

- the complete Feller-Pareto family and every Pareto II-IV alias
- inverse Burr, inverse paralogistic, loggamma, and several related aliases
- minimum-distance and empirical-moment estimation frameworks
- Hachemeister and general hierarchical credibility
- compound hierarchical simulation
- the full `ruin` fixed-point/Sparre-Andersen framework
- every vectorized/log-probability/tail flag supported by R's `d/p/q/r` API
- R S3 grouped-data classes, extraction methods, plotting, and histograms
- R conditions, formula handling, data frames, and serialization
- compiled use of the bundled dental and Hachemeister data objects
- integration through the upstream `expint` dependency

See `COVERAGE.md` for the detailed mapping and `PORTING_NOTES.md` for numerical
differences.

## License and provenance

The upstream package declares GNU GPL version 2 or later. This derivative work
uses `GPL-2.0-or-later` and includes the complete GPLv2 text in `LICENSE`.

- The unmodified upstream source is retained under `original/actuar-3.3-7`.
- The supplied archive is retained under `provenance`.
- SHA-256 manifests cover the supplied archive, retained source, and translated
  release files.
- Upstream authorship and contributor information is retained in the original
  `DESCRIPTION` file and summarized in `NOTICE`.

This project is an independent Fortran port and is not an official upstream
release.
