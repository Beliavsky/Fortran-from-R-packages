# Porting notes

## Translation strategy

`SpatialExtremes` mixes R orchestration with a sizeable C numerical backend. The Fortran port exposes reusable numerical models directly as array-oriented procedures rather than imitating R's S3/formula system.

The pairwise max-stable likelihoods follow the algebra of the upstream C routines. GEV margins are first transformed to unit Frechet and their log-Jacobians are added to each pair contribution, matching the package's composite-likelihood construction. Version 0.2.0 also exposes the contribution matrix itself, which makes score-based inference available without reconstructing the R object layer.

Exact finite-dimensional max-stable simulation is exposed separately from finite-storm spectral approximations. Gaussian-process conditioning is evaluated directly using the conditional mean and Schur-complement covariance; this is distributionally equivalent to the upstream unconditional-draw plus kriging-residual correction.

## Version 0.2.0 parity work

### Composite standard errors

`spatialextremes_standard_errors` implements the upstream Godambe/sandwich construction. A contribution callback is differentiated numerically, replicate scores form the variability matrix, pair-contribution scores form the sensitivity proxy used by upstream `standardErrors.c`, and the covariance is `H^{-1} J H^{-1}`. Stationary unit-Frechet wrappers are supplied for Smith, Schlather, Schlather-independence, Brown-Resnick, geometric-Gaussian, and extremal-t. Full GEV-margin contribution matrices are also public, allowing callers to include marginal/design parameters in the same generic calculation.

### Conditional max-linear simulation and partitions

The upstream `rcondMaxLin` hitting-scenario kernel is translated in `spatialextremes_maxlinear`. `spatialextremes_partitions` adds restricted-growth set-partition enumeration/canonicalization used by conditional max-stable simulation. Version 0.4.0 completes the numerical `condrmaxstab` path for Schlather, Brown-Resnick, and extremal-t: model-specific partition weights, exact enumeration for small conditioning sets, Gibbs updates for larger sets, upstream-style starting hitting scenarios, and extremal/sub-extremal conditional draws. The Brown-Resnick exact simulator also exposes its hitting scenario.

### Turning bands/random lines

`spatialextremes_turning_bands` translates the Van der Corput random-line construction, random rotations, and the 2D/3D turning-band Gaussian-process spectral kernels for Matern, Cauchy, powered exponential, Bessel, Gaussian, and fractional-Brownian models. Cartesian-grid mode is supported. Model-level Schlather, geometric-Gaussian, and extremal-t TBM simulators reproduce the upstream Poisson stopping construction and R-side default `uBound` values. Version 0.3.0 adds a self-contained radix-2 2-D FFT/circulant embedding engine and model-level Schlather, geometric-Gaussian, and extremal-t circulant simulators. The embedding uses powers of two rather than the upstream table of highly-composite FFT sizes, but implements the same covariance-embedding construction and automatically expands until the spectrum is nonnegative.

### Latent-variable kernels

The native DIC calculation, reusable GEV/Gaussian-field log-density pieces, and the full computational `latentgev` Metropolis-within-Gibbs chain driver are translated. The R formula/model-frame and S3 object orchestration remain intentionally outside the Fortran API.


## Version 0.4.0 parity closure

### Conditional max-stable workflow

`spatialextremes_mvprob` translates the randomized lattice integration used by the native conditional code for multivariate Normal and Student probabilities. `spatialextremes_conditional` uses these probabilities in Schlather, Brown-Resnick, and extremal-t block weights, exhaustive partition probabilities, and Gibbs partition updates. For larger conditioning sets, the starting partition follows the upstream strategy of taking a modal hitting scenario from unconditional simulations. The model-specific conditional simulators combine extremal functions associated with the selected partition and an independent sub-extremal Poisson process, and explicitly preserve the conditioning values.

The exact Brown-Resnick simulator was audited against upstream `rbrownexact`/`rhitscenbrown`; its log-scale maximum is initialized at the upstream negative sentinel rather than zero. This fixes an earlier porting error that imposed an unintended unit-Frechet lower bound of one.

### Design-matrix inference

`spatialextremes_design_inference` completes the standard-error convenience layer without duplicating the large analytic derivative blocks in `standardErrors.c`. Spatial location/scale/shape design matrices and optional observation-level temporal design matrices are expanded, transformed with the translated `gev2frechTrend` kernel, and differentiated through the existing per-pair contribution functions. `composite_sandwich_active` holds fixed parameters constant while differentiating only active entries. Generalized-Cauchy `smooth2` is included in the stationary and design-matrix standard-error parameter vectors.

## Helper reuse

The supplied MIT-licensed `r_mod.F90` is used for Normal/t probabilities and densities, RNGs, Bessel functions, and related R-compatible numerical helpers. BLAS/LAPACK are used for SPD Cholesky solves/inversion. No duplicate probability library was added.

## Optimizers

A compact Nelder-Mead routine is included for stationary likelihood and least-squares fits. Parameters whose domains are constrained are transformed internally so the optimizer does not repeatedly evaluate invalid covariance/shape values.

## Deliberate non-targets

Plotting, maps, S3 print/summary/profile methods, formula/model-frame parsing, datasets, and demo/presentation code are not translated. These are R interfaces rather than standalone numerical kernels.

## Validation

The retained tests cover the original numerical areas plus 0.2.0 regression tests for set partitions, conditional max-linear simulation, Van der Corput lines, point/grid TBM Gaussian simulation, TBM max-stable smoke tests, generic and Brown-Resnick sandwich standard errors, full-margin contribution decomposition, and latent DIC identities.
