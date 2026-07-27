# Computational coverage

This file maps the exported computational API of R package `sde` 2.0.21 to the
Fortran implementation. Names are modernized rather than copied verbatim.

## Exact model laws

| R API | Fortran API | Coverage |
|---|---|---|
| `rcOU`, `dcOU`, `pcOU`, `qcOU` | `ou_conditional_random`, `ou_conditional_pdf`, `ou_conditional_cdf`, `ou_conditional_quantile` | Complete scalar API |
| `rsOU`, `dsOU`, `psOU`, `qsOU` | `ou_stationary_random`, `ou_stationary_pdf`, `ou_stationary_cdf`, `ou_stationary_quantile` | Complete scalar API |
| `rcBS`, `dcBS`, `pcBS`, `qcBS` | `gbm_conditional_random`, `gbm_conditional_pdf`, `gbm_conditional_cdf`, `gbm_conditional_quantile` | Complete scalar API |
| `rcCIR`, `dcCIR`, `pcCIR`, `qcCIR` | `cir_conditional_random`, `cir_conditional_pdf`, `cir_conditional_cdf`, `cir_conditional_quantile` | Complete scalar API |
| `rsCIR`, `dsCIR`, `psCIR`, `qsCIR` | `cir_stationary_random`, `cir_stationary_pdf`, `cir_stationary_cdf`, `cir_stationary_quantile` | Complete scalar API |

Fortran elemental-style vectorization is obtained with ordinary loops or array
constructors. The scalar kernels avoid hidden allocation and can be called from
user vectorized wrappers.

## Simulation

| R API or `sde.sim` method | Fortran API | Coverage |
|---|---|---|
| `BM` | `brownian_motion` | Complete |
| `GBM` | `geometric_brownian_motion` | Complete |
| `BBridge` | `brownian_bridge` | Complete |
| `DBridge` | `diffusion_bridge_euler` | Complete numerical bridge kernel |
| Euler / predictor-corrector | `simulate_euler` | Complete, explicit `alpha`, `eta`, and derivative callback |
| Milstein | `simulate_milstein` | Complete |
| second Milstein | `simulate_milstein_second_order` | Complete |
| KPS | `simulate_kps` | Complete |
| Ozaki | `simulate_ozaki` | Complete |
| Shoji | `simulate_shoji` | Complete |
| conditional distribution | `simulate_conditional` | Complete through a sampler callback |
| exact OU | `simulate_ou_exact` | Complete |
| exact GBM | `simulate_gbm_exact` | Complete |
| exact CIR | `simulate_cir_exact` | Complete |
| exact acceptance method | `simulate_exact_ea` | Acceptance core complete; endpoint and `psi` supplied as callbacks |

## Transition densities and likelihoods

| R API | Fortran API | Coverage |
|---|---|---|
| `dcEuler` | `transition_density_euler` | Complete |
| `dcElerian` | `transition_density_elerian` | Complete |
| `dcKessler` | `transition_density_kessler` | Complete |
| `dcOzaki` | `transition_density_ozaki`, `ozaki_moments` | Complete |
| `dcShoji` | `transition_density_shoji`, `shoji_moments` | Complete |
| `dcSim` | `pedersen_transition_density` | Complete Pedersen Monte Carlo density |
| `EULERloglik` | `euler_log_likelihood` | Complete |
| `HPloglik` | `hermite_transition_density`, `hermite_log_likelihood` | Complete formula used by upstream C code |
| `SIMloglik` | `pedersen_log_likelihood` | Complete |
| exact model likelihoods | `ou_log_likelihood`, `gbm_log_likelihood`, `cir_log_likelihood` | Added convenience API |

## Estimation and testing

| R API | Fortran API | Coverage |
|---|---|---|
| `simple.ef` | `evaluate_simple_estimating`, `fit_simple_estimating` | Complete callback-based equivalent |
| `simple.ef2` | `evaluate_generator_estimating`, `fit_generator_estimating` | Complete callback-based equivalent |
| `linear.mart.ef` | `evaluate_linear_martingale`, `fit_linear_martingale` | Complete first/second-order callback equivalent |
| `gmm` | `gmm_moment_mean`, `gmm_hac_covariance`, `fit_gmm` | Complete two-stage iterative GMM structure |
| `sdeAIC` | `dc_transition_log_density`, `dc_log_likelihood`, `sde_aic`, `fit_sde_aic` | Complete numerical core |
| `sdeDiv` | `sde_divergence_test`, `fit_and_test_sde_divergence` | Complete test core and optional fit |

The R package delegates optimization to R's `optim`. The Fortran translation
uses an internal bounded Nelder-Mead implementation and exposes convergence
information in result types.

## Nonparametric and structural methods

| R API | Fortran API | Coverage |
|---|---|---|
| `ksdrift` | `kernel_drift` | Complete |
| `ksdiff` | `kernel_diffusion` | Complete |
| `ksdens` | `kernel_density` | Complete Gaussian KDE |
| `cpoint` | `detect_change_point` | Known-model and nonparametric modes |
| `MOdist` | `markov_operator_distance`, `bspline_basis_matrix` | Complete numerical construction; returns a symmetric matrix |

## Support routines added for a standalone package

The translation includes dependency-free normal, lognormal, gamma, chi-square,
and noncentral chi-square distributions; random variate generators; adaptive
quadrature; linear solves; matrix inversion; numerical Hessians; and bounded
optimization. These replace services supplied to the R package by R, MASS, and
other dependencies.

## Deliberate exclusions

The following are infrastructure rather than translated numerical algorithms:

- R expression parsing and symbolic differentiation.
- R `ts`, `zoo`, `fda`, and `dist` class construction.
- Plotting and print methods.
- Bundled `.rda` data sets and book scripts.
- R registration and `.Call` interface glue.

Users provide compiled callbacks, ordinary arrays, and explicit time steps in
the Fortran API.
