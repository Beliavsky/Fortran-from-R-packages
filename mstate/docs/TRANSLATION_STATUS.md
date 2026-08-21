# Translation status

Source: R package `mstate` 0.3.3 (2024-07-03).

## Implemented

| Upstream area | Fortran API | Status |
|---|---|---|
| `transMat`, `trans.comprisk`, `trans.illdeath` | `transition_from_matrix`, `trans_comprisk`, `trans_illdeath` | Computational core |
| `to.trans2`, `trans2Q`, `absorbing`, `is.circular` | `transition_map`, `trans2q`, `absorbing_states`, `is_circular` | Complete numerical equivalent |
| `trans2tra`, `tra2trans` | same names | Complete numeric conversion |
| `paths` | `enumerate_paths` | All reachable prefixes, matching upstream acyclic recursion |
| `msprep` / `msprepEngine` | `msprep` | Acyclic numeric matrix API; data-frame/formula plumbing omitted |
| `crprep` / `create.wData.omega` | `crprep` | Weighted expansion, censor/truncation KM weights, strata, keep matrix, shortening |
| `cutLMms` | `cut_landmark` | Core behavior |
| `xsect` | `xsect` | Core behavior |
| `events` | `event_counts` | Core counts/entering totals |
| `expand.covs.msdata` | `expand_covariates` | Numeric covariates; factor/model-matrix layer omitted |
| `msdata2etm`, `etm2msdata` | `msdata_to_etm`, `etm_to_msdata` | Core interval conversion |
| `msboot` | `bootstrap_msdata` | Subject-cluster resampling core |
| `agmssurv.c` | `agmssurv` | Direct numerical translation, Breslow and Efron |
| `msfit` post-Cox path | `msfit_from_cox_arrays` | Cox coefficients/covariance interface, offsets supported |
| `msfit` + `survival::coxph` | `coxph_fit_stratified_counting`, `msfit_from_survival_cox`, `msfit_cox` | Native stratified Cox bridge |
| `redrank` / `redrank.iter` | `redrank_fit` | Alternating stratified Cox fits and SVD normalization |
| `MarkovTest` | `markov_test` | Score traces, covariance, covariate adjustment, wild bootstrap, summaries |
| `probtrans` | `probtrans` | Forward/fixed-horizon product integral; Aalen/Greenwood covariance |
| `probtrans` with `Haz.boot` | `probtrans_bootstrap` | Point estimates plus bootstrap standard errors from hazard replicates |
| `ELOS` | `expected_length_of_stay` | Complete numerical integral |
| `LMAJ` | `landmark_aj` | Point estimates and single-/multi-start-state standard errors |
| `Cuminc` | `cumulative_incidence_fit`, `cumulative_incidence_grouped` | CIF/survival estimates, influence SEs, grouping |
| `mssample` / `mssample1` / `crsample` / `Hazsample` | `mssample_state`, `mssample_paths`, `mssample_data`, `sample_path_general` | Forward/reset clocks, history effects, censoring, all output modes |
| Markov-test weight helpers | `optimal_weights_multiple`, `optimal_weights_matrix` | Numerical formulas |
| `modify_transMat` | `modify_transition_relative` | Numeric state/transition splitting for absorbing targets |
| `varHaz.fixed` | `split_relative_hazards` | Population variance fixed at zero; excess inherits observed covariance |
| `haz_function` + `relsurv::expprep2` | `haz_function_relsurv` | Population/excess interval hazards, left truncation, event/add-time grid |
| `msfit.relsurv` | `msfit_relsurv`, `msfit_relsurv_full` | End-to-end rate-table hazard splitting; fixed/bootstrap/both workflows |
| `msboot.relsurv` / `.boot` | `msboot_relsurv` | Subject bootstrap, missing-transition handling, original/split hazards |
| `add.times` in `msfit.relsurv` | `hazard_add_times`, `add_times=` in `msfit_relsurv_full` | Carry-forward cumulative hazards/covariances |
| HLD/HMD rate-table loading | `relsurv_parsers::transrate_hld`, `transrate_hmd` | Included from relsurv-fortran dependency |

## Survival integration

The user-supplied `survival_f90.zip` is retained verbatim under
`vendor/survival_f90/`. The four modules needed by the mstate Cox path are also
compiled in `src/`: `survival_kinds`, `survival_types`, `survival_linalg`, and
`survival_cox`. Their LGPL-2.0-or-later SPDX declarations are retained.

The supplied survival implementation did not expose stratification in
`coxph_fit_counting`, so `mstate_cox` adds the separate-baseline/common-beta
stratified partial-likelihood layer required by `mstate` without modifying the
vendored survival sources.

## relsurv integration

`vendor/relsurv-fortran-v0.2.0.zip` retains the completed Fortran translation of
`relsurv` 2.3-3. `relsurv_kinds`, `relsurv_ratetable`, and `relsurv_parsers` are
compiled into this package. The relative-survival path therefore uses the
translated `pystep`/`pystep2` rate-table machinery and `expprep2_summary`
directly; population cumulative hazards no longer need to be supplied by the
caller.

The native API takes the rate-table variables as an already ordered numeric
matrix (`xrate`). R's `match.ratetable.mstate` and `rformulate.mstate` mostly
perform names, factor levels, date classes, formula parsing, model-frame
construction, and missing-value filtering. Those R runtime semantics are not
reproduced.

## Remaining parity boundaries

1. R missing-value/factor/formula/model-matrix behavior and S3 metadata.
2. R date/factor/name matching around ratetable objects. Native Fortran requires
   rate-table columns to be supplied in the table's dimension order.
3. Arbitrary R callback lists in `MarkovTest`; numerical traces are returned for
   native postprocessing instead.
4. Plotting/visualisation (`vis.mirror.pt`, `vis.multiple.pt`, plot/ggplot
   methods), intentionally excluded.
5. Native RNG streams are statistically equivalent but not bit-for-bit equal to
   R's RNG streams.

No substantial mstate-owned statistical algorithm remains intentionally
untranslated in this release.
