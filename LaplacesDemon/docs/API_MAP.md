# API map

## Optimization and approximation

Implemented numerical backends include:

- `bfgs_maximize`
- `newton_maximize` (NR)
- `levenberg_marquardt_maximize` (LM)
- `nelder_mead_maximize` (NM)
- `conjugate_gradient_maximize` (CG)
- `dfp_maximize`
- `lbfgs_maximize`
- `hooke_jeeves_maximize` (HJ)
- `hit_and_run_maximize` (HAR)
- `trust_region_maximize` (TR)
- `rprop_maximize`
- `sgd_maximize`
- `spg_maximize`
- `sr1_maximize`
- `pso_maximize`
- `genetic_maximize` (AGA)
- `soma_maximize`
- `bhhh_maximize`
- `laplace_approximation`
- `numerical_gradient`, `numerical_hessian`, `numerical_jacobian`
- `iterative_gauss_hermite`
- `componentwise_iterative_quadrature`
- `adaptive_sparse_grid_quadrature`
- `gauss_hermite_rule`, `gauss_hermite_cube`
- `variational_bayes_salimans2`

BHHH uses an explicit callback returning per-observation log-likelihood
components, replacing the implicit R model-list storage convention.

## MCMC

The upstream named algorithm catalog has fixed-dimensional numerical
counterparts:

| Upstream | Fortran API |
|---|---|
| ADMG | `admg_sample` |
| AFSS | `afss_sample` |
| AGG | `adaptive_griddy_gibbs_sample` |
| AHMC | `adaptive_hmc_sample` |
| AIES | `aies_sample` |
| AM | `adaptive_metropolis_sample` |
| AMM | `adaptive_mixture_metropolis_sample` |
| AMWG | `adaptive_mwg_sample` |
| CHARM | `charm_sample` |
| DEMC | `demc_sample` |
| DRAM | `dram_sample` |
| DRM | `drm_sample` |
| ESS | `elliptical_slice_sample` |
| GG | `griddy_gibbs_sample` |
| Gibbs | `gibbs_sample` |
| HARM | `hit_and_run_sample` |
| HMC | `hmc_sample` |
| HMCDA | `hmcda_sample` |
| IM | `independence_metropolis_sample` |
| INCA | `inca_sample` |
| MALA | `mala_sample` |
| MCMCMC | `mcmcmc_sample` |
| MTM | `multiple_try_metropolis_sample` |
| MWG | `mwg_sample` |
| NUTS | `nuts_sample` |
| OHSS | `ohss_sample` |
| pCN | `pcn_sample` |
| RAM | `ram_sample` |
| Refractive | `refractive_sample` |
| RDMH | `random_dive_sample` |
| RJ | `reversible_jump_selection_sample` |
| RSS | `rss_sample` |
| RWM | `rwm_sample` |
| SAMWG | `samwg_sample` |
| SGLD | `sgld_sample` |
| Slice | `slice_sample` |
| SMWG | `smwg_sample` |
| THMC | `tempered_hmc_sample` |
| twalk | `twalk_sample` |
| UESS | `uess_sample` |
| USAMWG | `usamwg_sample` |
| USMWG | `usmwg_sample` |

`Experimental` in upstream is intentionally not treated as a stable sampler.
INCA's cluster exchange is represented by an in-process multiple-chain
adaptation. Gibbs uses an explicit conditional-draw callback. RJ follows the
upstream fixed-vector variable-selection semantics with an active mask.

## Importance, Monte Carlo and evidence

- `importance_normal`
- `sir_normal`
- `rejection_normal`
- `pmc_normal`
- `bayesian_bootstrap_weights`
- `lml_laplace`
- `lml_harmonic`
- `lml_generalized_harmonic`

## Diagnostics

- acceptance rates
- IAT and ESS
- IMPS and batch-means MCSE
- Geweke
- Gelman-Rhat
- Raftery-Lewis
- Hangartner chi-square
- BMK
- Heidelberger-Welch
- KS
- WAIC and KLD

Plot-producing diagnostic wrappers are intentionally omitted; numerical
statistics are returned directly.

## Probability distributions and priors

The native library includes the scalar distributions from the earlier
releases plus the principal exported `distributions.R` families and
parameterization variants, including:

- Bernoulli/categorical/Dirichlet;
- Laplace, asymmetric/skew/log-Laplace and mixtures;
- half-normal/half-Cauchy/half-t;
- inverse-gamma/inverse-chi-square/inverse-Gaussian/inverse-beta;
- Pareto/generalized Pareto/generalized Poisson;
- scalar Student-t scale/precision CDF, quantile and RNG APIs;
- power exponential;
- multivariate normal, Student-t, Cauchy, Laplace and power-exponential in
  covariance, precision and Cholesky forms;
- Wishart/inverse-Wishart density/RNG/Cholesky forms;
- matrix normal, matrix gamma and inverse matrix gamma;
- normal-Laplace and normal/Laplace mixture priors;
- horseshoe, LASSO, Huang-Wand and Huang-Wand Cholesky;
- normal-Wishart and normal-inverse-Wishart;
- Yang-Berger and Cholesky form, hyper-g, Zellner;
- continuous-relaxation MRF;
- stick distribution;
- generic truncation density/CDF/quantile/RNG and truncated mean/variance.

The remaining unmatched `runifsphere` symbol in `distributions.R` is a local
nested helper used internally by RNGs, not an exported package API.

## Intentionally omitted R infrastructure

The following are not numerical omissions:

- plots and graphical diagnostics;
- print/summary methods and R classes/lists;
- formula/model-frame evaluation;
- monitor/output-file/database storage machinery;
- socket/cluster/HPC orchestration;
- dynamic package loading and R environment mutation.

Where those mechanisms carry an algorithmic callback, the Fortran API accepts
an explicit procedure and numeric arrays instead.
