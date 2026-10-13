# matrixStats API coverage

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

### Material per-family limitations

- Every operation accepts explicit real(dp) data. Convert integer and logical
  arrays with real(...,dp) and merge(1.0_dp,0.0_dp,...) respectively. Raw,
  character, factors, lists, R attributes and dimension-name handling are omitted.
- R NA and NaN are both represented by IEEE quiet NaN. Logical all/any results
  are real(dp): 1=true, 0=false, NaN=unknown, preserving three-valued logic.
- Rank `random` is omitted rather than using a different random generator.
- Set averages implement means (and a useful explicit weight-matrix extension),
  not arbitrary R FUN/tFUN callback dispatch. Each set is a column of integers;
  nonpositive entries are padding and are skipped.
- Most invalid numeric requests return NaN, missing-index sentinels, or empty
  arrays, rather than raising R errors/warnings. Shape contracts remain required.
- Scalar center broadcasting is not provided; matrix center arguments contain
  one center for every original row/column. R's periodic center validation is
  omitted. Use the correct sample mean/median when supplying a center.
- The singleton vector iqr and matrix IQRs deliberately differ on infinite
  inputs, matching upstream. Weighted mean, matrix weighted means, weighted
  variance/MAD and weighted median retain their different upstream weight rules.
- dp arithmetic is not bitwise-identical to R's extended-precision accumulators.
  Compensated/scaled sums and scaled direct products recover some extreme finite
  cases where upstream overflows; nonfinite mean residual corrections are skipped.

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

