# Computational coverage

Source package: `MarkowitzR` 1.0.2.0002.

## Exported routines

| R routine | Fortran interface | Coverage |
|---|---|---|
| `theta_vcov` | `theta_vcov` | Complete numerical translation |
| `itheta_vcov` | `itheta_vcov` | Complete numerical translation |
| `mp_vcov` | `mp_vcov` | Complete numerical translation |

## `theta_vcov`

Implemented behavior:

- fitted or omitted intercept
- lower-triangular packing of the unified second moment
- complete-case row deletion
- empirical covariance of the sample moment
- multivariate-normal analytic covariance
- native Bartlett/Newey-West HAC covariance
- custom covariance callback
- typed result with `mu`, `covariance`, `n`, and `pp`

The normal, no-intercept case is supported as an extension. Upstream R rejects
that combination because its specialized information-matrix path was only
written for the augmented moment.

## `itheta_vcov`

Implemented behavior:

- inversion of the estimated unified second moment
- lower-triangular packing of the inverse
- analytic delta method using the duplication matrix and Kronecker product
- empirical, normal, HAC, and callback covariance paths inherited from
  `theta_vcov`
- exact final symmetrization

## `mp_vcov`

Implemented behavior:

- unconditional portfolios
- conditional Markowitz coefficient matrices
- optional intercept
- optional observation weights
- upstream-compatible feature-only weighting
- explicit all-column weighting
- arbitrary subspace constraints through `jmat`
- arbitrary hedging constraints through `gmat`
- row-space validation for simultaneous `jmat` and `gmat`
- projected inverse second moments
- analytic delta-method covariance
- extraction of the coefficient matrix and its covariance

## Internal helpers

Translated or independently implemented equivalents include:

- `ivech`
- symmetric `vech`
- Kronecker products
- duplication matrices
- quadratic projection calculations
- projected precision matrices
- matrix inversion and rank calculation
- sample covariance and covariance of a sample mean

## Excluded infrastructure

The following are not computational algorithms and are not compiled:

- R package registration and namespace machinery
- R `lm` and `vcov` object dispatch
- optional `sandwich` package adapters
- roxygen/man pages, vignettes, and README plotting
- repository build, Docker, and bibliography scripts

Custom covariance estimation remains available through a typed Fortran callback.
