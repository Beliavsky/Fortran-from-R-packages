# Changelog

## 0.4.0

- Completed automatic conditional max-stable partition inference for Schlather, Brown-Resnick, and extremal-t, with exhaustive partition weights for small conditioning sets and Gibbs updates for larger sets.
- Added upstream-style starting hitting scenarios for all three conditional model families, including an exact Brown-Resnick hitting-scenario simulator.
- Added Brown-Resnick and extremal-t conditional extremal/sub-extremal simulation stages and high-level automatic conditional samplers.
- Added randomized-lattice/QMC multivariate Normal and Student CDF kernels used by the partition-weight calculations.
- Corrected Brown-Resnick exact simulation initialization so the log-scale maximum starts at the upstream negative sentinel rather than zero.
- Added `gev_to_frechet_trend` and a full spatial/temporal GEV design-matrix standard-error layer for all max-stable model families and spatial GEV.
- Added active/fixed-parameter sandwich calculations and generalized-Cauchy `smooth2` inference in stationary and design-matrix SE wrappers.
- Added `test_conditional_full.f90` and `test_design_inference.f90`; all ten test programs plus the example pass in the validation build.

## 0.3.0

- Added the full computational `latentgev` Metropolis-within-Gibbs driver: sitewise GEV updates, conjugate regression coefficients, inverse-gamma sills, lognormal range/smoothness MH updates, optional scale log-link, thinning/burn-in, chain storage, and acceptance/rejection diagnostics.
- Added a self-contained radix-2 two-dimensional FFT and circulant-embedding Gaussian-grid simulator. Embeddings grow automatically until the covariance spectrum is numerically nonnegative and are reused across storms.
- Added circulant max-stable simulation paths for Schlather, geometric-Gaussian and extremal-t processes.
- Added Schlather conditional max-stable simulation conditional on a supplied hitting partition, including conditional Student extremal-function simulation and the sub-extremal Poisson process.
- Added `test_remaining_gaps.f90`; all eight test programs pass with runtime checking and implicit-interface errors enabled.

## 0.2.0

- Added per-pair/per-replicate likelihood contribution matrices for Smith, Schlather, Schlather-independence, Brown-Resnick, geometric-Gaussian, and extremal-t, including full GEV margins.
- Added generic Godambe/sandwich standard-error machinery and stationary unit-Frechet model wrappers.
- Added conditional max-linear latent simulation and target-site conditional max-linear simulation.
- Added set-partition enumeration and canonicalization primitives for the conditional max-stable workflow.
- Added Van der Corput random lines, rotations, and 2D/3D turning-band Gaussian-process simulation with grid mode.
- Added Schlather, geometric-Gaussian, and extremal-t turning-band max-stable simulators.
- Added latent GEV likelihood, DIC, and Gaussian-field log-density kernels.
- Added `test_parity_targets` covering the new computational paths.
- Updated parity documentation to distinguish translated numerical kernels from the still-untranslated full latent MCMC, full `condrmaxstab`, and FFT/circulant paths.

## 0.1.0

- Initial modern Fortran/FPM computational port.
