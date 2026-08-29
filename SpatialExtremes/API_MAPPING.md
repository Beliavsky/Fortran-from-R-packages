# API mapping

This file maps the main computational R/native functionality to the modern Fortran modules. R-only formula, plotting, S3, and presentation routines are intentionally excluded.

| Upstream area | Fortran API/module |
|---|---|
| `dgev`, `pgev`, `qgev`, `rgev` | `spatialextremes_univariate::{dgev,pgev,qgev,rgev}` |
| `dgpd`, `pgpd`, `qgpd`, `rgpd` | `spatialextremes_univariate::{dgpd,pgpd,qgpd,rgpd}` |
| `gev2frech`, `gev2frechTrend` and inverse transforms | `gev_to_frechet`, `gev_to_frechet_trend`, `frechet_to_gev`, `gev_to_uniform` |
| `gevmle`, `gpdmle` numerical kernels | `spatialextremes_fit::{gevmle,gpdmle}` |
| covariance/correlation functions | `spatialextremes_covariance` |
| Smith Mahalanobis dependence | `mahalanobis_distances_2d/3d`, `extremal_coefficient_smith` |
| Schlather dependence | `extremal_coefficient_schlather`, `lplik_schlather` |
| Schlather + independent component | `extremal_coefficient_schlather_ind`, `lplik_schlather_ind` |
| Brown-Resnick | `brown_resnick_a`, `brown_resnick_loglik`, exact/spectral simulators |
| geometric Gaussian | `geom_gauss_a`, `geomgauss_loglik`, spectral/TBM simulators |
| extremal-t | `extremal_coefficient_extremalt`, `lplik_extremalt`, exact/spectral/TBM simulators |
| Smith pairwise likelihood | `lplik_smith`, `smith_loglik` |
| native pairwise contribution arrays | `lplik_*_contributions`, `*_loglik_contributions` |
| stationary `fitmaxstab` numerical fits | `fit_smith_frechet`, `fit_schlather_frechet`, `fit_brownresnick_frechet`, `fit_extremalt_frechet` |
| extremal-coefficient least squares | `spatialextremes_extcoeff_fit::lsfit_*` |
| `madogram`, `variogram`, `fmadogram`, concurrence kernels | `spatialextremes_dependence` |
| Gaussian process simulation | `simulate_gaussian_process` |
| `condrgp` numerical conditional Gaussian law | `conditional_gaussian_process` |
| `tbmcore`, `randomlines.c` | `simulate_gaussian_process_tbm`, `van_der_corput_lines`, `rotate_lines` |
| `rschlathertbm` | `simulate_schlather_tbm` |
| `rgeomtbm` | `simulate_geomgauss_tbm` |
| `rextremalttbm` | `simulate_extremalt_tbm` |
| other Schlather/geometric/extremal-t/Brown-Resnick simulation | `spatialextremes_simulation` |
| max-linear processes | `spatialextremes_maxlinear` |
| `rcondMaxLin` | `conditional_max_linear_latent`, `simulate_conditional_max_linear` |
| conditional max-stable set partitions | `list_set_partitions`, `canonicalize_partition`, `partition_block_count` |
| `condrmaxstab` partition weights | `schlather_partition_weights`, `extremalt_partition_weights`, `brownresnick_partition_weights` |
| `condrmaxstab` Gibbs partition updates | `gibbs_partitions_schlather`, `gibbs_partitions_extremalt`, `gibbs_partitions_brownresnick` |
| `condrmaxstab` starting hitting scenarios | `starting_partitions_schlather`, `starting_partitions_extremalt`, `starting_partitions_brownresnick` |
| conditional max-stable simulation | `sample_conditional_schlather`, `sample_conditional_extremalt`, `sample_conditional_brownresnick` and `*_given_partition` |
| native Gaussian/Student probability helpers used by conditional simulation | `mvnorm_cdf_qmc`, `mvstudent_cdf_qmc` |
| latent-model `DIC` native kernel | `latent_dic` |
| `latent` / native `latentgev` chain | `latent_gev_mcmc`, `latent_mcmc_result_t` |
| latent Gaussian-field prior kernel | `gaussian_field_logdensity` |
| Gaussian/Student copula likelihood and RNG | `spatialextremes_copula` |
| kriging core | `simple_kriging_weights`, `simple_kriging_predict` |
| spatial GEV likelihood core | `spatialextremes_spatgev` |
| `standardErrors.c` sandwich assembly | `composite_sandwich`, `composite_se_t` |
| stationary dependence SE wrappers | `*_frechet_standard_errors` |
| fixed/active parameter sandwich calculations | `composite_sandwich_active` and optional `active` masks |
| full spatial/temporal GEV design-matrix SE wrappers | `spatialextremes_design_inference::{smith_design_standard_errors,schlather_design_standard_errors,schlather_ind_design_standard_errors,brownresnick_design_standard_errors,geomgauss_design_standard_errors,extremalt_design_standard_errors,spatgev_design_standard_errors}` |
| circulant Gaussian process | `simulate_gaussian_grid_circulant` |
| Schlather `circ` simulation | `simulate_schlather_circulant` |
| geometric-Gaussian `circ` simulation | `simulate_geomgauss_circulant` |
| extremal-t `circ` simulation | `simulate_extremalt_circulant` |
| `AIC`, `TIC` numerical criteria | `aic_value`, `tic_value` |

## Remaining scope

No major standalone numerical target from the v0.3.0 parity list remains. The remaining omissions are R-side interface and presentation machinery: S3 classes, formula/model-frame construction, print/summary/plot/map/profile methods, `fields::Tps` spline presentation, and dataset/demo plumbing.

The Fortran implementation is not intended to be a byte-for-byte clone of R's internal implementation choices. In particular, circulant embedding uses a self-contained radix-2 FFT and the multivariate normal/Student probabilities used by conditional simulation use the translated randomized lattice/QMC algorithm.
