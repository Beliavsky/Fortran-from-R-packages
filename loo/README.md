# loo: Bayesian predictive assessment in modern Fortran

A substantial numerical translation of released **loo 2.10.1**, built with FPM.
The package provides Pareto-smoothed, truncated and standard importance sampling;
LOO/WAIC/ELPD; approximate-posterior correction; relative MCMC efficiency; weighted
predictive summaries; scoring; model comparison; stacking and pseudo-BMA; and fold
allocation. It is not an official Stan/loo release.

## Build and run

Place `loo/` at the root of Beliavsky/Fortran-from-R-packages, alongside the existing
`posterior/`, `matrixStats/`, and `rfortran-core/` translations. Then run:

```sh
cd loo
fpm build
fpm test
fpm run --example --all
```

A standard gfortran/FPM installation on Windows is sufficient. No separately
installed BLAS, LAPACK, R, or other system numerical library is required.
Sibling FPM path dependencies are reused without source copies. R is only needed
if explicitly regenerating the committed test-reference data. Tested with
gfortran 14.2.0 and FPM 0.13.0, default, checked-debug and release profiles
on Linux; a native Windows execution was not available.
No fast-math flags are supported: IEEE NaN/infinity behavior is part of the API.

## Main API and data layout

Import `loo_api` for the complete facade, or individual `loo_*` modules. All
maintained real arithmetic uses `dp` from `r_kinds`, re-exported by `loo_api`.
All numerical routines are pure except the callback entry points, which permit
caller-defined model evaluation. Every callback is a procedure pointer with
interface `log_likelihood_callback(observation, values)`; see the runnable example.
No global RNG or hidden mutable state is used.

```fortran
use loo_api
real(dp) :: log_lik(400, 8)
type(loo_result) :: assessment
type(importance_sampling_result) :: sampling
real(dp), allocatable :: weights(:, :)
! Fill log_lik with finite posterior-draw by observation log likelihoods.
assessment = loo(log_lik, r_eff=[1.0_dp])
sampling = psis(-log_lik)
weights = importance_weights(sampling, logarithm=.false.)
```

- Log-likelihood matrices are `(draw, observation)`. Three-dimensional arrays are
  `(iteration, chain, observation)` and flatten each chain contiguously, matching R.
- `psis`, `sis`, `tis` accept vectors, matrices and chain arrays. `r_eff` is an
  optional one-element array or one element per observation; NaN entries become 1.
  Positive finite efficiencies, including values above 1, are supported.
- `importance_sampling_result` contains raw and normalized log weights, log
  normalization constants, `pareto_k`, `n_eff`, `r_eff`, `tail_len`, method and status.
  Use `importance_weights`; do not reconstruct normalized weights by subtracting
  a very large raw log normalization constant.
- `loo_result%pointwise` has columns `(elpd, mcse_elpd, p_loo, looic, influence_k)`
  for LOO; `(elpd, p_waic, waic)` for WAIC; `(elpd, ic)` for ELPD. `estimates`
  has columns `(total, SE)`; LOO summarizes pointwise columns 1, 3 and 4.
- `relative_eff` requires **unlogged likelihoods**, not log likelihoods. Matrix and
  vector forms take equal-length chain IDs 1 through C, with interleaving supported.
  `relative_eff_function` takes an unlogged-likelihood callback.
- `E_loo` returns `expectation_result%value(probability_or_statistic, observation)`
  plus function-specific Pareto k. Quantiles require probabilities strictly inside
  (0,1). All-equal weights use R type 7; nonuniform weights interpolate the weighted CDF.
- `loo_compare` takes `(observation, model)` pointwise scores. Results are sorted
  best first, with original indexes in `order`, paired difference/SE, probability
  of worse performance, and each model's ELPD/SE. `compare` keeps the legacy
  two-model sign: second minus first.
- Predictive CRPS routines use the upstream **higher-is-better** sign convention.
  They accept explicit permutation columns; omitted permutations preserve pairing.
  LOO variants need ordinary normalized marginal and separately smoothed joint
  paired weights, not products of marginal smoothed weights.
- `stacking_weights` solves on the closed simplex and reports convergence status,
  objective and `kkt_gap`. `pseudobma_weights` without bootstrap rows is plain
  pseudo-BMA; supplied normalized `(bootstrap, observation)` rows produce
  pseudo-BMA+. `model_weights` routes these methods from supplied pointwise scores.
- Fold functions take caller-generated permutations and return one-based fold IDs.
  Their exact definitions, including the grouped mapping and R thinning precedence,
  are tested against the released R code. They do not generate random draws.

## Invalid inputs and numerical behavior

Derived results and model-weight routines report status 0 on success, 1 for invalid
arguments. Stacking status 2 means no certified convergence; predictive-score status
2 identifies a mathematically undefined scaled score or missing binary class.
Invalid results may leave arrays unallocated: check status before reading them.
`relative_eff` and invalid scalar numerical quantities return quiet NaN; invalid
fold allocations return zeros. Callback null pointers are rejected without calling.

All LOO/WAIC/ELPD likelihoods must be finite. WAIC needs at least two draws. Importance
log ratios may contain negative infinity (zero weight), but every column needs one
finite value; NaN and positive infinity are rejected. Too-short, constant or
unsuitable PSIS tails remain unsmoothed with k=positive infinity. A single observation
has an undefined between-observation standard error (NaN). Ordinary *raw*
weights can overflow when their mathematical exponential exceeds the real range.
Raw metadata at enormous offsets inevitably has real64 rounding; normalized weights
are stored independently in centered coordinates.

Documented numerical improvements over the release avoid loss of precision:
centered log-weight normalization; stable exponential tail differences; log-domain
MCSE that remains invariant to extreme likelihood shifts; centered weighted
variance; and overflow-safe tail selection for tiny relative efficiency. These
changes preserve the finite, well-conditioned reference outputs. See NOTICE.md.

## Translation coverage

Status: **substantial**. **40 of 44 (90.91%)** exported computational R functions
have meaningful mappings. This fraction measures mapped computational functions,
not complete R compatibility. There are 65 explicitly exported R names: 21
noncomputational/dispatch-only names are excluded with individual reasons in
API_COVERAGE.md. Registered but unexported S3 methods are not counted separately.
The following four computations remain untranslated:
`loo_moment_match`, `loo_moment_match.default`, `loo_subsample`, and
`loo_subsample.function`. They require complete model transformation/refitting or
subsampling workflows, beyond the supplied numerical arrays/callbacks here.

| R function | Fortran module | Fortran procedure(s) | Mapping status |
|---|---|---|---|
| `.compute_point_estimate` | `loo_evaluation` | `compute_point_estimate` | complete |
| `.thin_draws` | `loo_evaluation` | `thin_draws` | substantial |
| `E_loo` | `loo_diagnostics` | `E_loo` | substantial |
| `compare` | `loo_diagnostics` | `compare` | partial |
| `crps` | `loo_metrics` | `crps` | substantial |
| `elpd` | `loo_evaluation` | `elpd` | substantial |
| `kfold_split_grouped` | `loo_folds` | `kfold_split_grouped` | substantial |
| `kfold_split_random` | `loo_folds` | `kfold_split_random` | substantial |
| `kfold_split_stratified` | `loo_folds` | `kfold_split_stratified` | substantial |
| `loo` | `loo_evaluation` | `loo` | substantial |
| `loo.array` | `loo_evaluation` | `loo_array` | substantial |
| `loo.function` | `loo_evaluation` | `loo_function` | substantial |
| `loo.matrix` | `loo_evaluation` | `loo_matrix` | substantial |
| `loo_approximate_posterior` | `loo_evaluation` | `loo_approximate_posterior` | substantial |
| `loo_approximate_posterior.array` | `loo_evaluation` | `loo_approximate_posterior_array` | substantial |
| `loo_approximate_posterior.function` | `loo_evaluation` | `loo_approximate_posterior_function` | substantial |
| `loo_approximate_posterior.matrix` | `loo_evaluation` | `loo_approximate_posterior_matrix` | substantial |
| `loo_compare` | `loo_diagnostics` | `loo_compare` | substantial |
| `loo_crps` | `loo_metrics` | `loo_crps` | substantial |
| `loo_i` | `loo_evaluation` | `loo_i` | substantial |
| `loo_model_weights` | `loo_model_weights` | `model_weights` | partial |
| `loo_model_weights.default` | `loo_model_weights` | `model_weights` | partial |
| `loo_predictive_metric` | `loo_metrics` | `loo_predictive_metric` | substantial |
| `loo_scrps` | `loo_metrics` | `loo_scrps` | substantial |
| `mcse_loo` | `loo_diagnostics` | `mcse_loo` | complete |
| `pareto_k_ids` | `loo_diagnostics` | `pareto_k_ids` | complete |
| `pareto_k_table` | `loo_diagnostics` | `pareto_k_table` | substantial |
| `pseudobma_weights` | `loo_model_weights` | `pseudobma_weights` | substantial |
| `psis` | `loo_importance` | `psis` | substantial |
| `psislw` | `loo_importance` | `psislw` | substantial |
| `relative_eff` | `loo_evaluation` | `relative_eff`, `relative_eff_function` | substantial |
| `scrps` | `loo_metrics` | `scrps` | substantial |
| `sis` | `loo_importance` | `sis` | substantial |
| `stacking_weights` | `loo_model_weights` | `stacking_weights` | substantial |
| `tis` | `loo_importance` | `tis` | substantial |
| `waic` | `loo_evaluation` | `waic` | substantial |
| `waic.array` | `loo_evaluation` | `waic_array` | substantial |
| `waic.function` | `loo_evaluation` | `waic_function` | substantial |
| `waic.matrix` | `loo_evaluation` | `waic_matrix` | substantial |
| `weights.importance_sampling` | `loo_importance` | `importance_weights` | complete |

The manifest has one corresponding mapping entry per listed R function, with exact
source paths and material omissions. Plotting, printing, dynamic R classes, formula
processing, warning text, automatic parallel dispatch and implicit RNG are omitted.
The bare `kfold` generic is excluded because loo provides no computational method;
external modeling packages implement those methods.

## Tests and provenance

`fpm test` runs actual-R fixtures for importance weights, k/ESS, LOO, WAIC, ELPD,
approximate-posterior correction, E_loo (including diagnostics), deprecated psislw,
relative efficiency, comparisons, fold assignments, all five predictive metrics,
CRPS variants, stacking and pseudo-BMA. Independent tests cover chain ordering,
callbacks, NaN/Inf, tail ties, zero ratios, tiny efficiencies, constant samples,
large log offsets, variance offset invariance, classification and simplex boundaries.
Optional fixture scripts are in `test/reference/`; run them from this package root
with the pinned R package versions described in those scripts. Committed fixtures
let normal FPM tests run without R.

`example/predictive_assessment.f90` demonstrates end-to-end assessment, callbacks,
comparison, weighting and prediction. Full provenance, dependency decisions and
licenses are in NOTICE.md, LICENSE, LICENSE.upstream and upstream/SOURCE.json.
Original R sources are retained solely as provenance and are never compiled or
required by the Fortran runtime.
