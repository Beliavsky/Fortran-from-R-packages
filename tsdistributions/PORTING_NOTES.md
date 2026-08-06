# Porting notes

## Parameterization

The standardized location-scale parameterization follows upstream `tsdistributions`: `mu` is the mean and `sigma` is the standard deviation. Skew and shape meanings remain distribution-specific. GH and NIG conversion routines return `(mu, delta, beta, alpha)`.

## Numerical replacements

- TMB automatic differentiation and `nlminb` are replaced by bounded transformations, Nelder-Mead optimization, and finite-difference scores/Hessians.
- GH/NIG/GHST probability and random-generation support uses the compatible computational subset of the supplied `ghyp` Fortran translation.
- The SPD interior uses direct kernel sums rather than `KernSmooth::bkde`; tails use the upstream PWM GPD estimator.
- The authorized skewness-kurtosis boundary is obtained by deterministic constrained grid search rather than R interpolation and optimization helpers.
- Profiling is serial and deterministic for a fixed Fortran RNG seed; R `future` parallel infrastructure is omitted.

These choices preserve the statistical model but do not promise bit-for-bit agreement with TMB, R's RNG, `KernSmooth`, or external GH implementations.
