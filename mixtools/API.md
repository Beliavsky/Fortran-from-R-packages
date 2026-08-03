# API overview

All public interfaces are available from:

```fortran
use mixtools
```

R dots are written as underscores. For example, `regmixEM.lambda` becomes
`regmixEM_lambda`, and `test.equality.mixed` becomes
`test_equality_mixed`.

## Controls and results

`em_control` contains convergence tolerance, iteration limit, minimum scale,
ridge stabilization, seed, and verbosity settings.

Principal result types are:

- `mixture_result`
- `mv_mixture_result`
- `gamma_mixture_result`
- `multinomial_mixture_result`
- `regression_mixture_result`
- `semiparametric_result`
- `reliability_mixture_result`
- `model_selection_result`
- `bootstrap_result`
- `mcmc_result`

Each iterative result includes a status, convergence flag, iteration count,
log likelihood, and available log-likelihood history.

## Parametric mixture estimators

- `normalmixEM`
- `normalmixEM2comp`
- `normalmixMMlc`
- `tauequivnormalmixEM`
- `mvnormalmixEM`
- `gammamixEM`
- `multmixEM`
- `repnormmixEM`

Initialization-compatible one-step procedures are named `normalmix_init`,
`mvnormalmix_init`, `gammamix_init`, `multmix_init`, and
`repnormmix_init`.

## Regression mixtures

- `regmixEM`
- `regmixEM_lambda`
- `regmixEM_loc`
- `regmixEM_mixed`
- `logisregmixEM`
- `poisregmixEM`
- `segregmixEM`
- `hmeEM`
- `flaremixEM`
- `try_flare`
- `spregmix`
- `regmixMH`

The mixed-regression interface accepts integer group IDs and returns the
estimated group effects separately. The public implementation currently uses
ridge-regularized group intercepts.

## Semiparametric estimators

- `npEM`
- `npEMindrep`
- `npEMindrepbw`
- `npMSL`
- `spEM`
- `mvnpEM`
- `spEMsymloc`
- `spEMsymlocN01`

These return component densities on an explicit grid, posterior probabilities,
bandwidths, and locations where applicable.

## Reliability mixtures

- `expRMM_EM`
- `weibullRMM_SEM`
- `spRMM_SEM`

Inputs are event/censoring times and an integer event indicator (`1` for an
event, `0` for right censoring).

## Selection, bootstrap, and inference

- `normalmix_model_selection`
- `regmixmodel_sel`
- `multmixmodel_sel`
- `repnormmixmodel_sel`
- `boot_comp` / `normalmix_boot_comp`
- `boot_se` / `normalmix_boot_se`
- `test_equality`
- `test_equality_mixed`
- `post_beta`
- `regcr`

`test_equality` is generic for a univariate normal mixture or a regression
mixture. `boot_comp` and `boot_se` are currently the normal-mixture variants;
this is explicit in the typed API.

## Densities, simulation, and utilities

- `ddirichlet`, `dmvnorm`, `logdmvnorm`, `dexpmixt`
- `rnormmix`, `normmix_sim`, `rmvnorm`, `rmvnormmix`, `normmixrm_sim`
- `rexpmix`, `rweibullmix`, `rlnormscalemix`
- `wkde`, `wquantile`, `wIQR`
- `lambda`, `lambda_pert`
- `matsqrt`, `ellipse`, `depth`, `ldc`, `ldmult`
- `aug_x`, `makemultdata`, `perm`, `make_constraints`
- `compCDF`, `density_npEM`, `density_spEM`, `ise_npEM`
- `fdr_from_posterior`

The RNG is internal and seedable through `rng_state` and `rng_seed`, making
simulation and bootstrap runs deterministic for a given compiler-independent
integer arithmetic model.
