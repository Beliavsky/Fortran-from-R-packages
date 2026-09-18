# SteadyStateBVAR modern Fortran

This project is an unofficial modern free-form Fortran translation of the computational core of the R package **SteadyStateBVAR 0.2.0** by Mark Becker. It is intended to live as the top-level directory `SteadyStateBVAR/` in `Beliavsky/Fortran-from-R-packages`.

The translation implements the steady-state VAR parameterization, deterministic setup, OLS initialization, Minnesota and steady-state priors, beta restrictions, homoscedastic fitting, Random-Walk stochastic volatility, AR(1) stochastic volatility, joint predictive simulation, conditional homoscedastic forecasts, impulse responses, prior-probability-interval conversion, and posterior summaries. Plotting, S3/list presentation, `ts` metadata, console formatting, Stan objects, and other R-specific interface behavior are intentionally omitted.

## Build

Place this directory beside the shared repository packages `rfortran-core/`, `rfortran-linalg/`, and `MTS/`, then run:

```text
fpm build
fpm test
fpm run --example basic_workflow
fpm run --example stochastic_volatility
```

The manifest has no system BLAS/LAPACK link flags. It reuses:

- `rfortran-core` for the common `dp` kind, normal quantiles, and R type-7 quantiles;
- `rfortran-linalg` for checked dense solves, Cholesky factors, inverses, and SPD operations;
- the translated top-level `MTS` package for `var_psi_weights`, corresponding to `MTS::VARpsi`.

No dependency source is copied into this package.

## Public Fortran API

The umbrella module is `steadystatebvar`:

```fortran
use steadystatebvar, only : dp, bvar_model, forecast_result, irf_result
use steadystatebvar, only : bvar_summary_result
use steadystatebvar, only : bvar_create, setup_model, set_priors
use steadystatebvar, only : restrict_beta, fit_model, forecast_model
use steadystatebvar, only : conditional_forecast_model, irf_model
use steadystatebvar, only : ppi, bvar_summary
```

The numerical workflow mirrors the R package:

1. `bvar_create` copies an `(n,k)` data matrix into a `bvar_model`.
2. `setup_model` constructs lag/deterministic matrices and OLS quantities.
3. `set_priors` constructs Minnesota and steady-state priors plus homoscedastic or RW/AR1 stochastic-volatility priors.
4. `restrict_beta` applies the upstream near-zero convention by setting selected beta-prior variances to `1e-7_dp`.
5. `fit_model` samples the selected posterior and creates joint predictive draws.
6. `forecast_model`, `conditional_forecast_model`, and `irf_model` summarize predictive or response draws.
7. `bvar_summary` returns posterior mean or median parameter summaries; for SV fits its optional `t` uses the original-sample time index, as in R.
8. `ppi` converts a symmetric normal probability interval to implied mean and variance.

For stochastic volatility, pass `sv=.true.` and `sv_type="RW"` or `"AR1"` to `set_priors` together with the corresponding upstream prior arrays. `fit_model` accepts optional `sv_state_scale` and `sv_parameter_scale` proposal standard deviations. See `example/stochastic_volatility.f90`.

All maintained real variables use `dp` imported from `r_kinds`. Dummy arguments have explicit intent and FORD comments. The translation does not rely on `-ffast-math` or related finite-only assumptions.

## Bayesian fitting scope

Upstream `fit()` uses Stan/NUTS. The Fortran implementation targets the same model distributions but uses native MCMC:

- homoscedastic models use blocked Gibbs updates for `beta`, `Psi`, and `Sigma_u`;
- RW stochastic volatility uses Gaussian Gibbs updates for `beta`, `Psi`, and the free rows of `A`, scalar Metropolis updates for the latent log-volatility states, and inverse-gamma Gibbs updates for `phi`;
- AR(1) stochastic volatility uses Gaussian Gibbs updates for `beta`, `Psi`, and `A`, scalar Metropolis updates for latent log volatilities and `gamma_0`/`gamma_1`, and inverse-Wishart Gibbs updates for `Phi`.

The SV sampler uses the centered latent-state representation corresponding to the upstream Stan programs' noncentered parameterization. The likelihood, prior distributions, time-varying covariance construction, and predictive volatility recursion follow those Stan models. Posterior means/medians include the complete SV covariance path and the model-specific parameters (`A`, `phi`, or `gamma_0`, `gamma_1`, `Phi`).

This is not a Stan runtime translation. NUTS trajectories, adaptation, convergence diagnostics, chain objects, exact random-number streams, and Stan object construction are not reproduced. The SV Metropolis proposal scales are fixed rather than adaptively tuned, so practical mixing can differ materially from Stan even though the intended posterior target is the same.

The upstream R implementation explicitly rejects stochastic-volatility conditional forecasts. `conditional_forecast_model` does the same. For the supported homoscedastic case, the Fortran code uses the Gaussian projection form for equality-conditioned standard-normal shocks instead of the upstream SVD null-space construction; for a full-row-rank condition matrix the resulting conditional Gaussian distribution is equivalent.

## Translation coverage

Package status: **substantial**. **10 of 10 (100%)** exported computational R functions have a meaningful mapping. The fraction measures mapped computational functions, not complete R compatibility; Stan/NUTS execution details and R presentation layers remain outside the translation.

| R function | Fortran module | Fortran procedure(s) | Mapping status |
|---|---|---|---|
| `bvar` | `steadystatebvar_core` | `bvar_create` | substantial |
| `setup` | `steadystatebvar_core` | `setup_model` | substantial |
| `priors` | `steadystatebvar_core` | `set_priors` | substantial |
| `restrict_beta` | `steadystatebvar_core` | `restrict_beta` | substantial |
| `fit` | `steadystatebvar_core` | `fit_model` | substantial |
| `forecast` | `steadystatebvar_core` | `forecast_model` | substantial |
| `conditional_forecast` | `steadystatebvar_core` | `conditional_forecast_model` | substantial |
| `IRF` | `steadystatebvar_core` | `irf_model` | substantial |
| `ppi` | `steadystatebvar_core` | `ppi` | complete |
| `summary.bvar` | `steadystatebvar_core` | `bvar_summary` | substantial |

See `API_COVERAGE.md` for remaining differences and `PROVENANCE.md` for translation decisions.

## Validation

`test/test_steadystatebvar.f90` is deterministic and exercises model creation, deterministic specifications, priors, restrictions, homoscedastic Jeffreys/inverse-Wishart fits, RW-SV and AR1-SV fits, predictive simulation, conditional forecasts, time-specific OIRF/GIRF calculations, annualized growth transformations, and homoscedastic/SV summaries. The `example/` directory contains homoscedastic and RW-SV workflows.

See `VALIDATION.md` for the exact validation commands and the FPM availability note for the artifact-building environment.

## License and provenance

SteadyStateBVAR is GPL version 3 or later. This translation is therefore distributed under **GPL-3.0-or-later**. See `LICENSE`, `NOTICE.md`, `PROVENANCE.md`, and the retained upstream metadata/source files under `upstream/`.
