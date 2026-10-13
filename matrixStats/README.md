# matrixStats: modern Fortran numerical translation

A pure, free-form Fortran/FPM translation of matrixStats 1.5.0. Public procedures
are in `matrixstats`; the code reuses sibling `rfortran-core` for `r_kinds:dp`
and compatible sorting. No BLAS, LAPACK, R runtime, C compiler, network service,
or separately installed numerical library is needed to build or use it.

This is a modified independent translation, not an official matrixStats release.
License: Artistic-2.0 **and** GPL-2.0-or-later; see LICENSE, both full license texts,
NOTICE.md and UPSTREAM-NOTICES.txt. The small private R quicksort helper preserves
upstream weighted-median tie behavior; the ordinary shared stable sorter cannot
produce that behavior.

## Build, test and run

Place this directory at the root of Fortran-from-R-packages, alongside
`rfortran-core`, then run:

```text
cd matrixStats
fpm build
fpm test
fpm run --example
fpm build --profile release
fpm test --profile release
fpm run --example --profile release
```

The relative dependency is `rfortran-core = { path = "../rfortran-core" }`.
Only that sibling must be present; do not copy its sources here. The tested
repository revision is in NOTICE.md. Tested with FPM 0.13.0 and gfortran 14.2.0
on Linux, including debug bounds checking and release optimization. The sources
and manifest are Windows-compatible and use no system library links; native
Windows execution was not available for this validation.

Do **not** use `-ffast-math`, `-Ofast`, `-ffinite-math-only` or equivalent flags:
correct NaN/infinity semantics are part of this API.

## Basic use

```fortran
use matrixstats, only: dp, row_means2, row_quantiles, col_ranks
real(dp) :: x(2, 3)
real(dp), allocatable :: means(:), quartiles(:, :), ranks(:, :)
x = reshape([1.0_dp, 10.0_dp, 2.0_dp, 20.0_dp, 3.0_dp, 30.0_dp], [2, 3])
means = row_means2(x)
quartiles = row_quantiles(x, [0.25_dp, 0.5_dp, 0.75_dp])
ranks = col_ranks(x, ties_method = 'average', preserve_shape = .true.)
```

All maintained real data and constants use the single imported `dp` kind.
All computational procedures are pure; all array-valued operations allocate
results automatically. No mutable global state or random generator is used.
See `example/matrix_summary.f90` for a complete working example.

## API conventions

- Names use snake_case: `rowMeans2` becomes `row_means2`, `weightedMedian`
  becomes `weighted_median`, and `count` becomes `count_value` to avoid hiding
  the Fortran intrinsic. The complete mapping table appears below.
- Matrix data have normal Fortran/R column-major layout. Row reductions return
  one value per selected row; column reductions return one per selected column.
  Ranges and quantiles have summaries along their second dimension. `col_ranks`
  defaults to one original column per output row, matching R's preserveShape=FALSE;
  set `preserve_shape=.true.` for input orientation. Cumulative/difference
  outputs preserve the selected matrix orientation.
- Optional rows/cols/idxs use R numeric-index conventions: positive one-based
  entries select/reorder/repeat; negative entries exclude; zeros are discarded.
  Use `ms_na_index` for a missing index. Do not mix positive and negative indices.
  Missing or out-of-range positive selections yield NaN. Invalid mixed selectors
  yield a single missing selection instead of an R error.
- Use IEEE quiet NaN for either R NA or NaN. `na_rm` defaults false except
  `bin_means` (true) and set averages (automatic when missing data exist).
  Integer/logical data must be explicitly converted to real(dp). Names, storage
  modes, raw/character data, dimensions supplied separately from arrays, warning
  behavior, and R object construction are omitted.
- All/any predicates return real(dp) 0, 1 or NaN, so unknown truth is retained.
  Omitting `value` tests nonzero values, just like R's logical TRUE default;
  an explicit numeric `value` tests exact equality, and NaN tests missingness.
  Counts also use real(dp), allowing a NaN count. Missingness flags return logical.
- Quantiles support all nine Hyndman-Fan types (default 7). Ranks support
  max (default), min, average, first, last and dense. Random ties are omitted.
- Weighted medians interpolate by default; specifying `ties` changes that default
  to false, as in R. `interpolate=.true.` may explicitly override it. Ties accept
  weighted, min, max or mean. Zero/negative weights are dropped for median,
  variance and MAD; vector weightedMean follows upstream's actual signed-weight
  behavior while matrix weighted means and medians reject negative weights. See source dummy
  documentation and API_COVERAGE.md for missing/infinite-weight details.
- Matrix centers must have one entry per original row/column and are selected
  with the data. Scalar broadcasting and R's periodic center-validation warnings
  are omitted. A mismatched center length returns NaN summaries.
- `bin_means(y,x,bx,means,counts,...)` is a pure subroutine with allocated mean
  and optional count outputs. All other mapped operations are functions.
- `row_avgs_per_col_set` / `col_avgs_per_row_set` accept an integer matrix whose
  columns contain index sets, with nonpositive padding. They implement means,
  with optional same-shaped weight matrices; arbitrary R FUN/tFUN callbacks
  are omitted. Invalid weight shapes yield NaN summaries.
- Precondition violations generally yield NaNs, empty results or `ms_na_index`,
  not R errors. Data lengths and shapes must match documented dummy contracts.
  Indices use default Fortran integers, so arrays beyond that range are unsupported.
- Numerical results are tolerance-equivalent, not bitwise R clones. Compensated
  sums, overflow-rescaled sums/products and finite-only residual correction
  improve some extreme finite cases; they can differ from upstream overflow and
  accumulation order. NA-vs-NaN payload distinctions and signed-zero identity
  are not promised. No extra real kind or unsafe floating-point flag is used.

## Validation and regenerating references

The tests include deterministic R-generated references for every mapped export,
all nine quantile types, six deterministic tie methods, weighted tie policies,
subsets, repeated indices, missing/infinite data, boundary cases and shape
orientations. Additional tests check independent mathematical identities, empty
arrays, safe shape rejection and stable extreme arithmetic. No random seeds or
network access are needed when running `fpm test`.

To regenerate the oracle data, use R with **matrixStats 1.5.0** installed:

```text
python tools/generate_references.py
Rscript test/generate_order_reference.R
fpm test
```

The first command records auditable R expressions and embeds their reference
values in `test/test_reference.f90`; the second replaces `test/order_reference.txt`.
R is only needed for regenerating fixtures, never for normal build/test/use.
Run tests from the package root so the weighted oracle text file can be located.
See VALIDATION.md for exact commands and outcomes.

## Translation coverage

Package status: **substantial**. **99 of 99 (100%)** exported computational R
functions have meaningful mappings. This fraction measures mapped computational
functions, not complete R compatibility. There are no untranslated functions in
this basis. The upstream NAMESPACE has 102 exports: allocArray, allocMatrix and
allocVector are excluded because they are R runtime allocation/storage-interface
helpers, not computational algorithms. They are intentionally not implemented.
Alias exports such as rowAnyNAs are counted individually because they are
explicit exported computational entry points.

Rank functions and set-average functions are partial mappings for the reasons
below; all other mappings are substantial. No function claims complete R object,
error, storage-type, or interface compatibility.

| R function | Fortran module | Fortran procedure(s) | Mapping status |
|---|---|---|---|
| allValue | matrixstats | all_value | substantial |
| anyMissing | matrixstats | any_missing | substantial |
| anyValue | matrixstats | any_value | substantial |
| binCounts | matrixstats | bin_counts | substantial |
| binMeans | matrixstats | bin_means | substantial |
| colAlls | matrixstats | col_alls | substantial |
| colAnyMissings | matrixstats | col_any_missings | substantial |
| colAnyNAs | matrixstats | col_any_nas | substantial |
| colAnys | matrixstats | col_anys | substantial |
| colAvgsPerRowSet | matrixstats | col_avgs_per_row_set | partial |
| colCollapse | matrixstats | col_collapse | substantial |
| colCounts | matrixstats | col_counts | substantial |
| colCummaxs | matrixstats | col_cummaxs | substantial |
| colCummins | matrixstats | col_cummins | substantial |
| colCumprods | matrixstats | col_cumprods | substantial |
| colCumsums | matrixstats | col_cumsums | substantial |
| colDiffs | matrixstats | col_diffs | substantial |
| colIQRDiffs | matrixstats | col_iqr_diffs | substantial |
| colIQRs | matrixstats | col_iqrs | substantial |
| colLogSumExps | matrixstats | col_log_sum_exps | substantial |
| colMadDiffs | matrixstats | col_mad_diffs | substantial |
| colMads | matrixstats | col_mads | substantial |
| colMaxs | matrixstats | col_maxs | substantial |
| colMeans2 | matrixstats | col_means2 | substantial |
| colMedians | matrixstats | col_medians | substantial |
| colMins | matrixstats | col_mins | substantial |
| colOrderStats | matrixstats | col_order_stats | substantial |
| colProds | matrixstats | col_prods | substantial |
| colQuantiles | matrixstats | col_quantiles | substantial |
| colRanges | matrixstats | col_ranges | substantial |
| colRanks | matrixstats | col_ranks | partial |
| colSdDiffs | matrixstats | col_sd_diffs | substantial |
| colSds | matrixstats | col_sds | substantial |
| colSums2 | matrixstats | col_sums2 | substantial |
| colTabulates | matrixstats | col_tabulates | substantial |
| colVarDiffs | matrixstats | col_var_diffs | substantial |
| colVars | matrixstats | col_vars | substantial |
| colWeightedMads | matrixstats | col_weighted_mads | substantial |
| colWeightedMeans | matrixstats | col_weighted_means | substantial |
| colWeightedMedians | matrixstats | col_weighted_medians | substantial |
| colWeightedSds | matrixstats | col_weighted_sds | substantial |
| colWeightedVars | matrixstats | col_weighted_vars | substantial |
| count | matrixstats | count_value | substantial |
| diff2 | matrixstats | diff2 | substantial |
| indexByRow | matrixstats | index_by_row | substantial |
| iqr | matrixstats | iqr | substantial |
| iqrDiff | matrixstats | iqr_diff | substantial |
| logSumExp | matrixstats | log_sum_exp | substantial |
| madDiff | matrixstats | mad_diff | substantial |
| mean2 | matrixstats | mean2 | substantial |
| product | matrixstats | product | substantial |
| rowAlls | matrixstats | row_alls | substantial |
| rowAnyMissings | matrixstats | row_any_missings | substantial |
| rowAnyNAs | matrixstats | row_any_nas | substantial |
| rowAnys | matrixstats | row_anys | substantial |
| rowAvgsPerColSet | matrixstats | row_avgs_per_col_set | partial |
| rowCollapse | matrixstats | row_collapse | substantial |
| rowCounts | matrixstats | row_counts | substantial |
| rowCummaxs | matrixstats | row_cummaxs | substantial |
| rowCummins | matrixstats | row_cummins | substantial |
| rowCumprods | matrixstats | row_cumprods | substantial |
| rowCumsums | matrixstats | row_cumsums | substantial |
| rowDiffs | matrixstats | row_diffs | substantial |
| rowIQRDiffs | matrixstats | row_iqr_diffs | substantial |
| rowIQRs | matrixstats | row_iqrs | substantial |
| rowLogSumExps | matrixstats | row_log_sum_exps | substantial |
| rowMadDiffs | matrixstats | row_mad_diffs | substantial |
| rowMads | matrixstats | row_mads | substantial |
| rowMaxs | matrixstats | row_maxs | substantial |
| rowMeans2 | matrixstats | row_means2 | substantial |
| rowMedians | matrixstats | row_medians | substantial |
| rowMins | matrixstats | row_mins | substantial |
| rowOrderStats | matrixstats | row_order_stats | substantial |
| rowProds | matrixstats | row_prods | substantial |
| rowQuantiles | matrixstats | row_quantiles | substantial |
| rowRanges | matrixstats | row_ranges | substantial |
| rowRanks | matrixstats | row_ranks | partial |
| rowSdDiffs | matrixstats | row_sd_diffs | substantial |
| rowSds | matrixstats | row_sds | substantial |
| rowSums2 | matrixstats | row_sums2 | substantial |
| rowTabulates | matrixstats | row_tabulates | substantial |
| rowVarDiffs | matrixstats | row_var_diffs | substantial |
| rowVars | matrixstats | row_vars | substantial |
| rowWeightedMads | matrixstats | row_weighted_mads | substantial |
| rowWeightedMeans | matrixstats | row_weighted_means | substantial |
| rowWeightedMedians | matrixstats | row_weighted_medians | substantial |
| rowWeightedSds | matrixstats | row_weighted_sds | substantial |
| rowWeightedVars | matrixstats | row_weighted_vars | substantial |
| sdDiff | matrixstats | sd_diff | substantial |
| signTabulate | matrixstats | sign_tabulate | substantial |
| sum2 | matrixstats | sum2 | substantial |
| t_tx_OP_y | matrixstats | t_tx_op_y | substantial |
| varDiff | matrixstats | var_diff | substantial |
| weightedMad | matrixstats | weighted_mad | substantial |
| weightedMean | matrixstats | weighted_mean | substantial |
| weightedMedian | matrixstats | weighted_median | substantial |
| weightedSd | matrixstats | weighted_sd | substantial |
| weightedVar | matrixstats | weighted_var | substantial |
| x_OP_y | matrixstats | x_op_y | substantial |

## Explicit corrections of upstream defects

The following small, intentional differences are tested and disclosed rather
than silently reproducing upstream implementation defects:

1. With an odd number of infinite weights, upstream weightedMedian can read the
   wrong partially sorted element. For `x=[2,0,2], w=[Inf,Inf,Inf]`, matrixStats
   1.5.0 returns 0; this translation returns the mathematically correct median 2.
   The same repair applies to row/column weighted medians. Finite-weight tie
   ordering still follows R's exact unstable quicksort behavior.
2. Scalar weightedMad with no explicit weights ignores upstream `idxs` entirely;
   weightedVar/weightedSd with no explicit weights retain the wrong original
   weight-vector length after subsetting. Here all three consistently select
   observations first and then generate equal weights. For
   `x=[4,1,7,1,9,2], idxs=[1,2,6]`, this gives MAD=1.4826 and variance=7/3,
   versus upstream 2.9652 and 1.8666666667. Matrix wrappers are unaffected because
   their selection occurs before invoking vector summaries.
3. Compensated/scaled dp arithmetic recovers some finite results that upstream's
   specific accumulation order overflows. Mean residual corrections are skipped
   if they are nonfinite. These changes do not collapse actual NaNs or infinities.

