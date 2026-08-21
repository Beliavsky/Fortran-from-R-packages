# Parity notes for v0.4.0

## Relative-survival integration

Version 0.4.0 replaces the v0.3 external population-hazard matrix boundary with
the translated `relsurv` rate-table engine. The package compiles
`relsurv_kinds`, `relsurv_ratetable`, and `relsurv_parsers` from
`relsurv-fortran v0.2.0` and retains that full dependency archive under
`vendor/`.

`haz_function_relsurv` follows `mstate::haz_function`: rows are selected for one
transition; `Tstop`, event status, `Tstart`, and ordered population-table
covariates are passed to the translated `expprep2` fast path; the evaluation
grid is the union of follow-up times and requested additional times; and
interval hazards are

- population: `yidli / yi`
- excess: `dni / yi - yidli / yi`.

The reported standard error follows the source cumulative
`sqrt(cumsum(dni / yi^2))` expression. `population_hazards_relsurv` then
accumulates population increments onto the original `msfit` time grid.
`msfit_relsurv` deliberately computes cumulative excess hazard as observed
cumulative hazard minus cumulative population hazard, matching the source's more
numerically stable route.

The Fortran `xrate` matrix corresponds to the numeric `R` matrix produced after
R's `rformulate.mstate`/`match.ratetable.mstate` processing. Column ordering is
explicit. Date classes, factor-name matching, and formula evaluation are R
interface features and are not emulated.

`time_format='days'|'years'|'months'` in `msfit_relsurv_full` uses the upstream
365.241-day year and 365.241/12-day month conversions. Rate-table covariates
such as age/calendar-date are assumed already expressed in the units expected
by the supplied rate table, as they are after R's `rmap` evaluation.

## Relative-survival bootstrap

`msboot_relsurv` follows `msboot.relsurv`/`msboot.relsurv.boot` numerically:
subjects are sampled with replacement; the no-covariate stratified Cox fit is
represented by the equivalent per-transition Breslow/Nelson-Aalen estimator;
replicate hazards are carried onto the requested original time grid; and the
relative-survival split is recomputed through the rate table for every
replicate.

A transition absent from one bootstrap sample is omitted for that transition,
not treated as a failure of the whole replicate. Pointwise sample variances use
all available replicates for that transition. Replicates that contain only one
distinct source subject are excluded, matching the source safeguard. Optional
`boot_original=.true.` also retains variances for the original unsplit hazards,
corresponding to `boot_orig_msfit`.

`probtrans_bootstrap` implements the `Haz.boot` branch of upstream `probtrans`:
each valid cumulative-hazard replicate is propagated through the same product
integral and the standard deviation at each time/start/end-state combination is
returned as the transition-probability SE. Missing transitions in a replicate
have zero hazard, which is equivalent to the upstream bootstrap object omitting
that transition.

R RNG streams are not reproduced bit-for-bit. The resampling mechanism and
sample-variance formulas are the same.

## Cox bridge

`survival_f90.zip` supplies a counting-process Cox fitter but no `strata`
argument. `mstate` requires stratified Cox models with a common regression
coefficient vector and separate baseline hazards. `coxph_fit_stratified_counting`
therefore sums the Cox score and information over strata while preserving the
Breslow/Efron tie formulas used by the supplied survival code.

`msfit_cox` groups rows by stratum before calling the translated `agmssurv`
kernel. `msfit_from_survival_cox` is also available when a fitted
`coxph_result` already exists.

## `redrank`

The Fortran API accepts numeric design matrices rather than R formulas. It
reproduces the two alternating Cox regressions from `redrank.iter`, the
transition-specific score expansion, and the final SVD normalization of
`Alpha` and `Gamma` while leaving `Beta = Alpha*Gamma` unchanged.

The upstream R source accepts a `method` argument in `redrank`, but the current
`redrank.iter` calls shown in mstate 0.3.3 do not pass that argument to `coxph`.
The Fortran API intentionally honors its `method` argument; use
`method='efron'` to reproduce the effective default Cox tie method of that R
source path.

## `MarkovTest`

The numerical implementation follows `R/MarkovTest.R`, including the
Cox-adjusted score term, covariance calculation and wild bootstrap. The
covariate-adjusted covariance expression preserves the source's use of the
first qualifying-state `Zbar` index in the second factor where that appears in
the R code, rather than silently correcting it.

Wild-bootstrap random streams are native Fortran streams. Distributional
semantics are centered Poisson(1) or standard normal multipliers.

## `crprep`

The implementation follows `crprep.R` and `create.wData.omega.R` at the numeric
array level. Event-time grids and product-limit weights are constructed
separately within each stratum. Competing-event subjects remain in the expanded
risk set at later event times, while subjects with the event of interest stop
at their own event. The small machine-precision shifts used upstream to order
tied failures, censorings, and entries are retained.

## `Cuminc` and `LMAJ`

`cumulative_incidence_fit` computes Aalen-Johansen competing-risks estimates
and standard errors from subject case-weight influence derivatives. Grouped
curves are fitted independently by distinct numeric group.

For multiple landmark starting states, `landmark_aj` includes uncertainty in
both empirical start-state frequencies and within-start-state AJ estimates.

## `paths` and `mssample`

`enumerate_paths` returns every reachable path prefix, including the starting
state by itself and nonterminal prefixes, matching upstream `paths()`.

The `mssample` translation supports forward and reset clocks, visit-history
state times, state-history regression effects, censoring, and state/path/data
output. Native RNG streams are not bit-identical to R.

## Bounds-safety correction found during v0.4 testing

The relative-survival bootstrap exercised a previously untested ordering of
event times and exposed a Fortran short-circuit assumption in the v0.3
`sort_unique` insertion loop. It used a compound `j>=1 .and. y(j)>key`
condition. Fortran does not guarantee short-circuit evaluation, so the loop was
rewritten with an explicit boundary test. This changes no statistical formula
but makes the Nelson-Aalen path bounds-clean under `-fcheck=all`.
