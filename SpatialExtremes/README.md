# SpatialExtremes-fortran

Modern Fortran/FPM computational port of the R package **SpatialExtremes 2.1-0**.

The port focuses on reusable numerical/statistical kernels rather than R's S3, formula, plotting, and map interfaces. It uses standard free-form Fortran and links to BLAS/LAPACK. The supplied `r_mod.F90` is reused for R-compatible probability, random-number, and special-function helpers.

## Included computational areas

- GEV and GPD density/CDF/quantile/RNG, log likelihoods, and GEV <-> unit-Frechet/uniform transforms.
- Distance utilities and covariance/correlation models used by SpatialExtremes: powered exponential, Whittle-Matern, Cauchy/generalized Cauchy, Bessel, fractional-Brownian covariance, and anisotropic Smith Mahalanobis distances.
- Smith, Schlather, Schlather-with-independence, Brown-Resnick/geometric-Gaussian, and extremal-t dependence kernels.
- Pairwise composite log likelihoods for Smith, Schlather, Schlather-independence, Brown-Resnick/geometric-Gaussian, and extremal-t models, including optional pair weights and GEV marginal Jacobians.
- Per-replicate/per-pair contribution matrices for both unit-Frechet and full GEV-margin likelihoods. These support composite score and sandwich calculations without an R object layer.
- Univariate GEV/GPD MLE and stationary max-stable pairwise-likelihood fits for Smith, Schlather, Brown-Resnick, and extremal-t models.
- Empirical madogram, variogram, F-madogram/extremal-coefficient transforms, concurrence estimators, lambda-madograms, and Smith/Schlather-Tawn empirical extremal-coefficient estimators.
- Least-squares extremal-coefficient fitting for Smith, Schlather, Schlather-independence, Brown-Resnick, geometric-Gaussian, and extremal-t models.
- Gaussian-process simulation and conditioning/kriging.
- Spectral simulation for Schlather, geometric-Gaussian, extremal-t, and Brown-Resnick models.
- Exact finite-dimensional simulation for Schlather, extremal-t, and Brown-Resnick models.
- Turning-band/random-line Gaussian-process simulation in 2D/3D, including Cartesian-grid mode, plus Schlather, geometric-Gaussian, and extremal-t TBM max-stable simulators.
- Radix-2 FFT/circulant-embedding Gaussian-grid simulation plus Schlather, geometric-Gaussian, and extremal-t circulant max-stable simulators.
- Max-linear evaluation, unconditional simulation, Smith-style grid design matrices, and conditional latent/max-linear simulation following the upstream hitting-scenario kernel.
- Full conditional max-stable numerical workflow for Schlather, Brown-Resnick, and extremal-t: set partitions, model-specific partition weights, exhaustive sampling for small conditioning sets, Gibbs partition updates for larger sets, upstream-style starting partitions/hitting scenarios, and conditional extremal/sub-extremal simulation.
- Full latent-GEV Metropolis-within-Gibbs chain driver with conjugate beta/sill updates and range/smoothness Metropolis steps.
- Gaussian and Student copula likelihood/simulation kernels.
- Spatial-GEV design-matrix parameter expansion and likelihood kernels.
- Generic composite-likelihood sandwich standard errors, stationary unit-Frechet wrappers, active/fixed-parameter support, and full spatial/temporal GEV design-matrix standard-error wrappers for Smith, Schlather, Schlather-independence, Brown-Resnick, geometric-Gaussian, extremal-t, and spatial GEV models.
- Latent-model GEV likelihood, DIC, and Gaussian-field log-density kernels.
- AIC/TIC numerical helpers, finite-difference Hessian, Stirling/Bell/Vandermonde and penalization helpers.

The public umbrella module is `SpatialExtremes`.

## Build

With FPM installed:

```text
fpm test
fpm run --example demo
```

The project links to BLAS and LAPACK.

## Scope boundary

Version 0.4.0 closes the two specialized numerical parity targets documented in 0.3.0. The `condrmaxstab` computational workflow is now represented for Schlather, Brown-Resnick, and extremal-t, including automatic hitting-partition inference and model-specific conditional simulation. The standard-error layer now accepts spatial and temporal GEV design matrices and an active-parameter mask corresponding to upstream fixed-parameter handling. Generalized-Cauchy `smooth2` is included as a differentiated parameter where applicable.

The intentionally omitted material is R interface/presentation code rather than standalone numerical kernels: S3/formula/model-frame construction, print/summary/plot/map/profile methods, dataset/demo plumbing, and `fields::Tps` presentation. The Fortran FFT uses an independent radix-2 implementation instead of translating R's bundled `fft.c`, while preserving the circulant-embedding computation.

See `API_MAPPING.md`, `PORTING_NOTES.md`, and `CHANGELOG.md` for details.

## License

SpatialExtremes-derived code is GPL-2.0-or-later, following the upstream package. `r_mod.F90` remains under its supplied MIT license. The upstream source and copyright material used for provenance is retained under `upstream/`.
