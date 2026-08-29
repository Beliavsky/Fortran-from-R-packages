# API and computational coverage

## Modern public API

The `kfas` module exposes a typed interface around the main numerical paths:

- `type(kfas_model)` -- numerical state-space arrays and exact-diffuse metadata.
- `type(kfas_filter_result)` -- predicted and filtered states/covariances,
  innovations, gains, optional original-scale signal, and log-likelihood.
- `type(kfas_smooth_result)` -- smoothed states, signals, observation and state
  disturbances, associated variances, and smoothing recursions.
- `type(kfas_approx_result)` -- Gaussian approximation, linear predictor,
  pseudo-observations, pseudo-variances, convergence information, and status.
- `kfas_gaussian_filter` -- ordinary and exact-diffuse Gaussian filtering.
  Correlated multivariate observation covariance `H` is automatically converted
  to the sequential-processing form using the upstream LDL transformation.
- `kfas_gaussian_smooth` -- Gaussian state, signal, and disturbance smoothing.
- `kfas_gaussian_loglik` -- Gaussian log-likelihood through the same filtering
  path, including automatic correlated-`H` handling.
- `kfas_ldl_transform` -- explicit KFAS LDL observation-equation transformation.
- `kfas_approximate_nongaussian` -- Gaussian approximation for Gaussian,
  Poisson, binomial, gamma, and negative-binomial observation families.
- `kfas_nongaussian_loglik` -- deterministic Laplace/approximation likelihood
  (`nsim = 0` in the R API terminology). The importance-sampling kernels are
  retained in the library but are not hidden behind a random-number policy.
- `kfas_init_theta` -- native version of KFAS `initTheta`.
- Distribution codes `kfas_gaussian`, `kfas_poisson`, `kfas_binomial`,
  `kfas_gamma`, and `kfas_negative_binomial`.
- `kfas_ar_transform` -- KFAS AR-parameter recursion.
- `kfas_ldl_factor` -- KFAS LDL decomposition.
- `kfas_weighted_mean_cov` -- normalized weighted mean/covariance calculation.

For an LDL-transformed correlated `H`, smoothed state and signal quantities are
reported on the original state/signal scale. `obs_disturbance` follows the
internally transformed sequential observation equation; the result flag
`observation_disturbance_transformed` records this explicitly.

## Retained computational kernels

All 30 upstream `src/*.f90` numerical source units are included in the FPM
library and preserve the corresponding KFAS algorithms:

| Area | Retained source units |
| --- | --- |
| Gaussian filtering / likelihood | `filter1step`, `filter1stepnovar`, `kfilter`, `kfilter2`, `gloglik` |
| Gaussian smoothing | `gsmoothall`, `smoothonestep` |
| Non-Gaussian approximation | `approx`, `approxloop`, `ngfilter`, `ngloglik`, `ngsmooth`, `ptheta`, `pytheta` |
| Simulation / importance sampling | `filtersimfast`, `simfilter`, `simgaussian`, `simgaussianuncond`, `smoothsim`, `smoothsimfast`, `isample`, `isamplefilter` |
| Signal/state helpers | `kfstheta`, `predict`, `mvfilter`, `marginalxx` |
| Linear algebra / transformations | `ldl`, `ldlssm`, `artransform`, `covmeanw` |

The specialized simulation, importance-sampling, non-Gaussian filtering, and
non-Gaussian smoothing routines remain compiled with their upstream numerical
calling conventions. They are retained rather than replaced by reduced or
approximate reimplementations.

## R-layer code intentionally not reproduced as R APIs

`SSModel`, `SSMtrend`, `SSMseasonal`, `SSMcycle`, `SSMarima`, `SSMregression`,
`SSMcustom`, and `SSMbespoke` mainly build R list/formula objects and metadata.
Fortran callers populate `kfas_model` arrays directly instead of recreating
those R object conventions.

`fitSSM` delegates generic optimization to R's `optim`; an unrelated optimizer
is not embedded in this package. Fortran callers can optimize
`kfas_gaussian_loglik` or `kfas_nongaussian_loglik` with their preferred
optimizer.

S3 methods for printing, plotting, extraction, time-series attributes, formula
processing, and presentation-oriented residual handling are intentionally
omitted. Plotting and interactive functionality are explicitly out of scope.
