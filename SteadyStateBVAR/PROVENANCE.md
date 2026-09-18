# Provenance

## Source material

The translation was prepared from the attached source archive `SteadyStateBVAR-master.zip`, whose package metadata identifies SteadyStateBVAR version 0.2.0 and Mark Becker as author, maintainer, and copyright holder. Retained upstream material is under `upstream/`.

Computational exports were identified from the upstream `NAMESPACE`. Plotting-only exports `steady_state_priors_plot` and `stochastic_volatility_plot` were intentionally excluded from the coverage denominator under the repository translation convention.

The numerical translation used these upstream files directly, retained as repository-relative provenance copies:

- `upstream/R/bvar.R`
- `upstream/R/setup.R`
- `upstream/R/priors.R`
- `upstream/R/restrict_beta.R`
- `upstream/R/fit.R`
- `upstream/R/forecast.R`
- `upstream/R/conditional_forecast.R`
- `upstream/R/irf.R`
- `upstream/R/ppi.R`
- `upstream/R/summary.R`
- `upstream/stan/steady_state_bvar_homoscedastic_jeffreys_prior.stan`
- `upstream/stan/steady_state_bvar_homoscedastic_inverse_wishart_prior.stan`
- `upstream/stan/steady_state_bvar_RW_stochastic_volatility.stan`
- `upstream/stan/steady_state_bvar_AR1_stochastic_volatility.stan`
- `upstream/stan/include/license.stan`

The `upstream/stan/` files are copies of the corresponding upstream package `inst/stan/` files; the shorter retained path avoids reproducing generated RStan C++ sources.

## Shared dependency decision

Before implementing reusable numerical helpers, the target repository was checked for compatible top-level translations. Existing `rfortran-core`, `rfortran-linalg`, and `MTS` projects supply the needed kind/distribution/quantile, dense linear-algebra, and VAR psi-weight APIs, respectively. The translation therefore references those sibling projects in `fpm.toml` and contains no copied dependency source.

`rfortran-linalg` owns the provenance and license obligations of its LAPACK backend. SteadyStateBVAR does not call LAPACK directly and does not vendor or declare its own BLAS/LAPACK implementation.

## Homoscedastic posterior translation

The homoscedastic model is implemented as a blocked Gibbs sampler. It targets the same Gaussian steady-state VAR likelihood and the same independent normal priors on `vec(beta)` and `vec(Psi)` used by the upstream Stan models. Conditional on `beta` and `Psi`, `Sigma_u` is sampled from the inverse-Wishart conditional implied by either the upstream Jeffreys prior or the upstream inverse-Wishart prior.

## Stochastic-volatility posterior translation

Both upstream SV specifications are now represented numerically.

For RW stochastic volatility:

- `A` is the same unit-lower-triangular contemporaneous matrix, with free parameters ordered row-major as upstream;
- `log_lambda(1,:)` uses the upstream independent normal initial-state prior based on the diagonal of `Omega_log_lambda_1`;
- latent log variances follow driftless random walks;
- each `phi_i` has the upstream inverse-gamma prior and is sampled from its centered full conditional;
- time-varying `Sigma_u,t` and predictive `Sigma_u,pred` use `A^{-1} Lambda_t A^{-T}` exactly as in the Stan generated quantities.

For AR(1) stochastic volatility:

- `gamma_0` and `gamma_1` use the upstream independent normal priors based on covariance diagonals, with `gamma_1` restricted to `(-1,1)`;
- correlated log-volatility innovations use `Phi` with the upstream inverse-Wishart prior;
- latent log variances follow `log_lambda_t = gamma_0 + gamma_1 * log_lambda_{t-1} + nu_t`;
- in-sample and predictive covariance matrices use the same construction as upstream.

The Stan programs use a noncentered latent `z` parameterization. The Fortran sampler works in the corresponding centered `log_lambda` coordinates. The resulting centered transition densities include the parameter-dependent normalizing terms through the inverse-gamma/inverse-Wishart full conditionals, so the intended posterior distribution is the coordinate-transformed counterpart of the upstream model. The transition implementation is nevertheless a different MCMC algorithm and is not claimed to reproduce Stan draws or finite-run behavior.

Homoscedastic beta/Psi/Sigma blocks, SV beta/Psi/A blocks, RW `phi`, and AR1 `Phi` use Gibbs draws. Latent SV states and AR1 `gamma_0`/`gamma_1` use symmetric random-walk Metropolis updates. Proposal scales are fixed caller inputs rather than Stan-style adaptive tuning.

## Forecast, IRF, and summary decisions

Predictive recursions follow the upstream Stan generated quantities, including future RW or AR1 log-volatility evolution. The fit object retains raw time-varying covariance draws plus posterior mean/median covariance paths, matching the numerical content of upstream `fit()` summaries.

For `conditional_forecast`, the upstream SVD null-space draw is replaced by the mathematically equivalent full-row-rank Gaussian projection

`eta = R' (R R')^-1 r + [I - R' (R R')^-1 R] z`,

where `z` is standard normal. The upstream implementation explicitly refuses stochastic-volatility conditional forecasts, and the Fortran procedure retains that behavior.

For SV `IRF` and `bvar_summary`, optional `t` is interpreted as an original-sample time index. It is converted internally to the estimated-sample covariance index by subtracting `p`, matching the R convention. Default `t` is the final original observation.

## Deliberate non-parity

This is a model/computational translation, not a Stan runtime or R presentation translation. NUTS trajectories, Stan adaptation and diagnostics, multi-chain Stan objects, exact R/Stan RNG streams, S3/list construction, names, `ts` metadata, printing, rounding, parameter-filter presentation, and plotting are intentionally not reproduced.
