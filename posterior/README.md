# posterior: modern Fortran computational translation

FPM translation of **posterior 1.7.0** for numeric posterior draws, convergence
diagnostics, generalized Pareto fitting and predictive checks. License:
BSD-3-Clause. See [NOTICE.md](NOTICE.md) for pinned provenance and authors.

## Build and use

Place this directory at the repository root alongside `rfortran-core`:

```text
Fortran-from-R-packages/
  posterior/
  rfortran-core/
```

Then run from `posterior/`:

```sh
fpm build
fpm test
fpm run --example
```

The only runtime dependency is the sibling FPM path package `rfortran-core`.
No installed R, BLAS, LAPACK or ARPACK library is needed. The package uses
standard free-form Fortran and a single `dp` kind imported from `r_kinds`.
It is designed for gfortran/FPM, including Windows; the recorded test platform
is Linux, so no native Windows execution claim is made.

Import everything with `use posterior`, or import a focused module:

- `posterior_diagnostics`: rank/folded split Rhat, nested Rhat, basic/bulk/tail/
  mean/SD/quantile ESS, MCSE, ranks, quantiles, autocovariance and nine default
  draw summaries
- `posterior_pareto`: generalized-Pareto CDF/quantiles and weighted fitting,
  tail smoothing, log-weight smoothing, k diagnostics and sample-size helpers
- `posterior_predictive`: weighted empirical and Pareto-refined PIT,
  POT/PRIT/PIET uniformity tests, entropy, dissention and numeric modes
- `posterior_draws`: numeric array conversion, binding, subsetting, chain
  splitting/merging, thinning, reordering, resampling, log-weight columns,
  within/across-draw reductions, CDFs, conditional broadcasting and matrix products
- `posterior_transforms`: numeric rvar matching/casting, matrix diagonals,
  callback-based mutation, per-draw actions and random-matrix margin reductions

```fortran
use posterior, only: dp, rhat, ess_bulk, ess_tail, mcse_mean
real(dp) :: x(1000, 4)  ! Populate with iterations by chains.
! ... fill x ...
print *, rhat(x), ess_bulk(x), ess_tail(x), mcse_mean(x)
```

See `example/diagnose_draws.f90` for a complete deterministic program. Every
procedure's dummy arguments have FORD documentation in the source, including
array shape, optional behavior and status meaning. Public numerical procedures
are pure or pure elemental, except `for_each_draw`, which permits the caller's
explicit side-effecting action.

## Array and error conventions

- A scalar variable across chains is `x(iteration, chain)`. A draws collection
  is `x(iteration, chain, variable)`. Flattening runs iterations first, then chains.
- A numeric rvar vector is `x(draw, element)`; a random matrix is
  `x(draw, row, column)`. Array values carry no hidden names or R class metadata.
- Diagnostics return IEEE NaN when undefined. Public subroutines with `status`
  return zero for success and documented nonzero values for invalid or
  unsuitable data. Never compile with fast-math, finite-math-only or analogous
  assumptions: NaN/infinity behavior is part of the API.
- Callback transforms require an associated procedure pointer implementing
  the published interface. `mutate_variables` appends/replaces one scalar
  variable per call; `rvar_apply_margin` retains one matrix margin and reduces
  each draw separately. They do not evaluate R text or formulas.
- `summarise_draws` returns rows in this order: mean, median, SD, MAD, q5, q95,
  Rhat, bulk ESS, tail ESS. The `summarize_draws` R alias maps to that procedure.
- `rvar_summary` and `rvar_over_draws_summary` select an operation using its
  documented method string. Their axis directions intentionally differ.
- Pareto routines accept flattened scalar draws. Their automatic efficiency
  treats these as one chain. For chain-aware `pareto_smooth`, `pareto_khat` and
  `pareto_diags`, compute `ess_tail(x)/size(x)` on the chain matrix and pass
  `r_eff`. Distribution parameters must be finite with positive scale.
- `pit` and `pareto_pit` take one predictive-variable column at a time.
  Explicit `uniform` controls randomization within ties; its default is 0.5
  (mid-PIT), unlike R's automatic random draw. Resampling similarly accepts
  explicit uniforms for reproducibility. No attempt is made to reproduce R RNG.
- `weight_draws_array` appends or replaces an explicit log-weight column;
  `extract_draws_weights` retrieves/normalizes it. Pareto-smooth weights
  separately when desired. No automatic R attribute or variable-name discovery.
- `quantile2` implements type 7. It preserves exact infinite endpoints;
  missing-value removal is optional. Numeric missing values all use IEEE NaN,
  without R's separate `NA` payload and logical/integer missing types.
- R's pinned diagnostic split behavior is preserved: odd chains lose their
  center; with only 2 or 3 iterations, R's dropped single-row dimensions produce
  a chains-by-two matrix. The public `split_draws_chains` instead discards the
  last odd iteration, matching the exported R operation. These are different
  upstream operations, not interchangeable shortcuts.
- The upstream Geyer ESS cap is preserved without warning text. Biased
  autocovariance is O(N²) here, via shared direct computation rather than FFT.
- The pinned `ps_convergence_rate` has rate zero at k=0 and k=1; this behavior is
  retained. Invalid inputs are explicitly guarded rather than relying on R errors.
- `is_constant` is a logical finite-range test with no missing logical state;
  Nonfinite/empty inputs and invalid tolerances return false; very large ranges
  are checked without overflow. Numeric category summaries cannot represent
  unobserved factor levels. These limitations are reflected in mapping notes.

## Verification

`fpm test` includes deterministic numerical properties and checked-in reference
fixtures generated by actual **R 4.5.0 / posterior 1.7.0**. The fixture scripts
are reproducible from the package root:

```sh
Rscript test/reference_diagnostics.R
Rscript test/reference_pareto.R
Rscript test/reference_draws.R
Rscript test/reference_transforms.R
```

R is needed only to regenerate/recheck references. Ordinary Fortran tests run
without R. Tests include autocorrelated and antithetic chains, shifted/scaled
chains, ties, odd and short chains, single and degenerate chains, missing and
infinite draws, weighted/zero-weight GPD fits, log weights, PIT tails, matrix
broadcasting, missing reductions, callbacks and exact distribution properties.
See [VALIDATION.md](VALIDATION.md) for the final recorded commands and results.

## Translation coverage

Package status: **substantial**. **91 of 113 (80.53%)** exported computational R
functions have meaningful mappings. The fraction measures mapped computational
functions, including partial and substantial specializations, rather than
complete R compatibility. Exactly 27 of the 140 explicit exports are excluded
as metadata/class/name/default-list or example-data interfaces. Untranslated
numeric conversions, constructors and transformations remain in the denominator.

The 22 untranslated computational exports are `as_draws`, `as_draws_df`,
`as_draws_list`, `as_draws_rvars`, `as_rvar_factor`, `as_rvar_numeric`,
`as_rvar_ordered`, `draws_df`, `draws_list`, `draws_matrix`, `draws_rvars`, `drop`,
`extract_list_of_variable_arrays`, `extract_variable_array`, `rdo`,
`repair_draws`, `rfun`, `rstar`, `rvar`, `rvar_factor`, `rvar_ordered`, `rvar_rng`.

[API_COVERAGE.md](API_COVERAGE.md) records every inclusion/exclusion.
`fpm.toml` contains exact upstream paths and material limitations per mapping.
Registered but not explicitly exported S3 method names are not counted again.

| R function | Fortran module | Fortran procedure(s) | Mapping status |
|---|---|---|---|
| %**% | posterior_draws | rvar_matmul | substantial |
| %in% | posterior_transforms | rvar_in | substantial |
| E | posterior_draws | rvar_expectation | substantial |
| Pr | posterior_draws | rvar_probability | substantial |
| as_draws_array | posterior_draws | matrix_to_draws | substantial |
| as_draws_matrix | posterior_draws | draws_to_matrix | substantial |
| as_rvar | posterior_draws | rvar_constant | partial |
| as_rvar_integer | posterior_transforms | rvar_as_integer | substantial |
| as_rvar_logical | posterior_transforms | rvar_as_logical | substantial |
| autocorrelation | posterior_diagnostics | autocorrelation | substantial |
| autocovariance | posterior_diagnostics | autocovariance | substantial |
| bind_draws | posterior_draws | bind_draws_arrays | substantial |
| cdf | posterior_draws | rvar_cdf | substantial |
| diag | posterior_transforms | rvar_diagonal, rvar_diagonal_matrix | substantial |
| dissent | posterior_predictive | dissent | substantial |
| draws_array | posterior_draws | matrix_to_draws | substantial |
| entropy | posterior_predictive | entropy | substantial |
| ess_basic | posterior_diagnostics | ess_basic | substantial |
| ess_bulk | posterior_diagnostics | ess_bulk | substantial |
| ess_mean | posterior_diagnostics | ess_mean | substantial |
| ess_median | posterior_diagnostics | ess_median | substantial |
| ess_quantile | posterior_diagnostics | ess_quantile | substantial |
| ess_sd | posterior_diagnostics | ess_sd | substantial |
| ess_tail | posterior_diagnostics | ess_tail | substantial |
| extract_variable | posterior_draws | extract_draws_variable | substantial |
| extract_variable_matrix | posterior_draws | extract_draws_variable_matrix | substantial |
| for_each_draw | posterior_transforms | for_each_draw | partial |
| gpdfit | posterior_pareto | gpdfit | substantial |
| is_constant | posterior_diagnostics | is_constant | partial |
| mad | posterior_draws | rvar_over_draws_summary | partial |
| match | posterior_transforms | rvar_match | substantial |
| mcse_mean | posterior_diagnostics | mcse_mean | substantial |
| mcse_median | posterior_diagnostics | mcse_median | substantial |
| mcse_quantile | posterior_diagnostics | mcse_quantile | substantial |
| mcse_sd | posterior_diagnostics | mcse_sd | substantial |
| merge_chains | posterior_draws | merge_draws_chains | substantial |
| modal_category | posterior_predictive | modal_category | substantial |
| mutate_variables | posterior_transforms | mutate_variables | partial |
| order_draws | posterior_draws | order_draws_array | substantial |
| pareto_convergence_rate | posterior_pareto | pareto_convergence_rate | substantial |
| pareto_diags | posterior_pareto | pareto_diags | substantial |
| pareto_khat | posterior_pareto | pareto_khat | substantial |
| pareto_khat_threshold | posterior_pareto | pareto_khat_threshold | substantial |
| pareto_min_ss | posterior_pareto | pareto_min_ss | substantial |
| pareto_pit | posterior_predictive | pareto_pit | substantial |
| pareto_smooth | posterior_pareto | pareto_smooth | substantial |
| pgeneralized_pareto | posterior_pareto | pgeneralized_pareto | substantial |
| pit | posterior_predictive | pit | substantial |
| ps_convergence_rate | posterior_pareto | ps_convergence_rate | substantial |
| ps_khat_threshold | posterior_pareto | ps_khat_threshold | substantial |
| ps_min_ss | posterior_pareto | ps_min_ss | substantial |
| ps_tail | posterior_pareto | ps_tail | substantial |
| ps_tail_length | posterior_pareto | ps_tail_length | substantial |
| qgeneralized_pareto | posterior_pareto | qgeneralized_pareto | substantial |
| quantile2 | posterior_diagnostics | quantile2 | substantial |
| r_scale | posterior_diagnostics | r_scale | substantial |
| resample_draws | posterior_draws | resample_draws_array, subset_draws_array | partial |
| rhat | posterior_diagnostics | rhat | substantial |
| rhat_basic | posterior_diagnostics | rhat_basic | substantial |
| rhat_nested | posterior_diagnostics | rhat_nested | substantial |
| rvar_all | posterior_draws | rvar_summary | substantial |
| rvar_any | posterior_draws | rvar_summary | substantial |
| rvar_apply | posterior_transforms | rvar_apply_margin | partial |
| rvar_ifelse | posterior_draws | rvar_ifelse | substantial |
| rvar_is_finite | posterior_draws | rvar_special_masks | substantial |
| rvar_is_infinite | posterior_draws | rvar_special_masks | substantial |
| rvar_is_na | posterior_draws | rvar_special_masks | partial |
| rvar_is_nan | posterior_draws | rvar_special_masks | partial |
| rvar_mad | posterior_draws | rvar_summary | substantial |
| rvar_max | posterior_draws | rvar_summary | substantial |
| rvar_mean | posterior_draws | rvar_summary | substantial |
| rvar_median | posterior_draws | rvar_summary | substantial |
| rvar_min | posterior_draws | rvar_summary | substantial |
| rvar_prod | posterior_draws | rvar_summary | substantial |
| rvar_quantile | posterior_draws | rvar_quantiles | substantial |
| rvar_range | posterior_draws | rvar_ranges | substantial |
| rvar_sd | posterior_draws | rvar_summary | substantial |
| rvar_sum | posterior_draws | rvar_summary | substantial |
| rvar_var | posterior_draws | rvar_summary | substantial |
| sd | posterior_draws | rvar_over_draws_summary | partial |
| split_chains | posterior_draws | split_draws_chains | substantial |
| subset_draws | posterior_draws | subset_draws_array | substantial |
| summarise_draws | posterior_diagnostics | summarise_draws | partial |
| summarize_draws | posterior_diagnostics | summarise_draws | partial |
| thin_draws | posterior_draws | thin_draws_array | partial |
| u_scale | posterior_diagnostics | u_scale | substantial |
| uniformity_test | posterior_predictive | uniformity_test | substantial |
| var | posterior_draws | rvar_over_draws_summary | partial |
| variance | posterior_draws | rvar_over_draws_summary | partial |
| weight_draws | posterior_draws | weight_draws_array, extract_draws_weights | partial |
| z_scale | posterior_diagnostics | z_scale | substantial |
