# API coverage

Coverage basis: distinct exported computational R functions. Plotting,
printing, formatting, datasets, the general `dlm()` R-list/S3 object constructor,
S3 getters/setters, `is.dlm`/`as.dlm`, `dropFirst` time-series metadata handling,
and registered-but-unexported methods are not in the denominator.

Package status: **substantial**. Mapped functions: **25 / 25 = 100%**.

## Main numerical coverage

- Model constructors: regression, polynomial trend, seasonal, trigonometric,
  and univariate ARMA.
- Model composition: binary additive latent model and binary outer sum.
- State-space inference: Gaussian filtering, missing observation components,
  time-varying system matrices, constant-free negative log likelihood,
  RTS smoothing, one-step residuals, and deterministic forecasting.
- Sampling: backward Gaussian state paths, the `dlmGibbsDIG` precision Gibbs
  workflow, Wishart simulation, and random constant DLM generation.
- Convex support geometry and sampling: `convex.bounds` plus the Gilks ARMS
  envelope/Metropolis algorithm, including the upstream multivariate
  random-direction hit-and-run wrapper.
- Maximum likelihood: `dlm_mle` builds candidate models through a typed callback
  and minimizes `dlm_ll` through the sibling `optimx` package.
- Utilities: stationary AR coefficient transformation, covariance reconstruction,
  Sokal MCMC standard errors, means, and ergodic cumulative means.

## Material differences

- `dlmFilter`/`dlmSmooth` return full covariance arrays instead of R lists of
  SVD factors. The implemented algebra is covariance-form Kalman/RTS rather
  than a line-by-line copy of the upstream square-root SVD kernels.
- `dlmModARMA` currently covers the univariate branch only.
- `dlmForecast` returns forecast moments; R's `sampleNew` branch is omitted.
- `dlmRandom` is a meaningful constant-model generator but does not reproduce
  the upstream time-varying generator or exact transition-matrix draw law.
- Random-number streams are Fortran intrinsic RNG streams and are not expected
  to match R bit-for-bit. ARMS therefore preserves the algorithm rather than the
  exact R random sequence.
- `dlm_mle` replaces R's arbitrary `build` closure, `...`, and `optim` result list
  with typed builder/result objects. Its optimizer methods are those provided by
  the repository's top-level `optimx` translation. The callback bridge is
  synchronous and intentionally rejects reentrant `dlm_mle` calls.
- `bdiag` is represented by a typed two-matrix `block_diag2` routine rather than
  an arbitrary heterogeneous R list/variadic interface.

## Coverage completion

There are no unmapped functions in the stated exported-computational-function
coverage basis. A 100% mapped-function fraction does not imply complete R/S3
interface compatibility: several mappings remain substantial or partial for the
reasons above and in `fpm.toml`.

The detailed one-entry-per-mapped-function metadata is in `fpm.toml`.
