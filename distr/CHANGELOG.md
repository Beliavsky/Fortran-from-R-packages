# Changelog

## 0.2.0

- Added dependency-free radix-2 FFT utilities.
- Added `convolve_fft` for continuous grid convolution and lattice mass convolution.
- Added `convpow_fft`; `convpow(...,grid_points=...)` now selects the FFT path when appropriate.
- Added closed-form convolution simplifications for additional compatible families.
- Added closed-form affine simplifications for several named families.
- Added `sf`, `logcdf`, and `logsf`, plus vector entry points.
- Added direct upper-tail quantiles, including logarithmic tail probabilities.
- Added stable upper-tail formulas for named central and noncentral laws.
- Added direct log-density evaluation for major named distributions and mixtures.
- Added `log_transform`, `sqrt_transform`, `reciprocal_transform`, and `power_transform`.
- Added `lattice_dist` and `weighted_empirical_dist` constructors.
- Improved left-limit CDF behavior for mixtures and composed distributions with atoms.
- Added raw/central moments, skewness, and excess kurtosis methods.
- Added `test_v02` regression coverage for the new numerical paths.

## 0.1.0

- Initial standalone Fortran/FPM computational translation of `distr` 2.9.7.
