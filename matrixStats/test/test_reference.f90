! SPDX-License-Identifier: Artistic-2.0
! Generated deterministic references: R 4.5.0, matrixStats 1.5.0.
program test_reference
   use matrixstats
   use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_positive_inf, ieee_is_nan
   implicit none
   real(dp) :: x(4, 5), xn(4, 5), v(6), vn(6), w(6), wr(5), wc(4)
   real(dp) :: nan, inf, ninf
   real(dp), allocatable :: bm(:)
   integer :: checks, failures
   checks = 0
   failures = 0
   nan = ieee_value(0.0_dp, ieee_quiet_nan)
   inf = ieee_value(0.0_dp, ieee_positive_inf)
   ninf = -inf
   x = reshape([1.0_dp, 4.0_dp, -2.0_dp, 9.0_dp, 2.0_dp, 4.0_dp, -1.0_dp, 8.0_dp, &
      3.0_dp, 8.0_dp, 0.0_dp, 7.0_dp, 4.0_dp, 8.0_dp, 1.0_dp, 6.0_dp, 5.0_dp, 9.0_dp, 2.0_dp, 5.0_dp], [4, 5])
   xn = x
   xn(1, 2) = nan
   xn(3, 4) = nan
   v = [4.0_dp, 1.0_dp, 7.0_dp, 1.0_dp, 9.0_dp, 2.0_dp]
   vn = v
   vn(2) = nan
   w = [1.0_dp, 2.0_dp, 0.0_dp, 3.0_dp, 4.0_dp, 2.0_dp]
   wr = [1.0_dp, 0.0_dp, 2.0_dp, 3.0_dp, 1.0_dp]
   wc = [1.0_dp, 2.0_dp, 0.0_dp, 3.0_dp]
   call bin_means(vn, v, [0.0_dp, 2.0_dp, 4.0_dp, 8.0_dp, 10.0_dp], bm)
   call check('rowSums2', row_sums2(x), &
      [real(dp) :: 15.0_dp, 33.0_dp, 0.0_dp, 35.0_dp])
   call check('rowSums2 selected', row_sums2(xn, rows = [4, 1, 4], cols = [5, 2, 3], na_rm = .true.), &
      [real(dp) :: 20.0_dp, 8.0_dp, 20.0_dp])
   call check('rowSums2 missing', row_sums2(xn), &
      [real(dp) :: nan, 33.0_dp, nan, 35.0_dp])
   call check('rowMeans2', row_means2(x), &
      [real(dp) :: 3.0_dp, 6.5999999999999996_dp, 0.0_dp, 7.0_dp])
   call check('rowMeans2 selected', row_means2(xn, rows = [4, 1, 4], cols = [5, 2, 3], na_rm = .true.), &
      [real(dp) :: 6.666666666666667_dp, 4.0_dp, 6.666666666666667_dp])
   call check('rowMeans2 missing', row_means2(xn), &
      [real(dp) :: nan, 6.5999999999999996_dp, nan, 7.0_dp])
   call check('rowVars', row_vars(x), &
      [real(dp) :: 2.5_dp, 5.7999999999999998_dp, 2.5_dp, 2.5_dp])
   call check('rowVars selected', row_vars(xn, rows = [4, 1, 4], cols = [5, 2, 3], na_rm = .true.), &
      [real(dp) :: 2.333333333333333_dp, 2.0_dp, 2.333333333333333_dp])
   call check('rowVars missing', row_vars(xn), &
      [real(dp) :: nan, 5.7999999999999998_dp, nan, 2.5_dp])
   call check('rowSds', row_sds(x), &
      [real(dp) :: 1.5811388300841898_dp, 2.4083189157584592_dp, 1.5811388300841898_dp, 1.5811388300841898_dp])
   call check('rowSds selected', row_sds(xn, rows = [4, 1, 4], cols = [5, 2, 3], na_rm = .true.), &
      [real(dp) :: 1.5275252316519465_dp, 1.4142135623730951_dp, 1.5275252316519465_dp])
   call check('rowSds missing', row_sds(xn), &
      [real(dp) :: nan, 2.4083189157584592_dp, nan, 1.5811388300841898_dp])
   call check('rowMedians', row_medians(x), &
      [real(dp) :: 3.0_dp, 8.0_dp, 0.0_dp, 7.0_dp])
   call check('rowMedians selected', row_medians(xn, rows = [4, 1, 4], cols = [5, 2, 3], na_rm = .true.), &
      [real(dp) :: 7.0_dp, 4.0_dp, 7.0_dp])
   call check('rowMedians missing', row_medians(xn), &
      [real(dp) :: nan, 8.0_dp, nan, 7.0_dp])
   call check('rowMads', row_mads(x), &
      [real(dp) :: 1.4825999999999999_dp, 1.4825999999999999_dp, 1.4825999999999999_dp, 1.4825999999999999_dp])
   call check('rowMads selected', row_mads(xn, rows = [4, 1, 4], cols = [5, 2, 3], na_rm = .true.), &
      [real(dp) :: 1.4825999999999999_dp, 1.4825999999999999_dp, 1.4825999999999999_dp])
   call check('rowMads missing', row_mads(xn), &
      [real(dp) :: nan, 1.4825999999999999_dp, nan, 1.4825999999999999_dp])
   call check('rowIQRs', row_iqrs(x), &
      [real(dp) :: 2.0_dp, 4.0_dp, 2.0_dp, 2.0_dp])
   call check('rowIQRs selected', row_iqrs(xn, rows = [4, 1, 4], cols = [5, 2, 3], na_rm = .true.), &
      [real(dp) :: 1.5_dp, 1.0_dp, 1.5_dp])
   call check('rowIQRs missing', row_iqrs(xn), &
      [real(dp) :: nan, 4.0_dp, nan, 2.0_dp])
   call check('rowMins', row_mins(x), &
      [real(dp) :: 1.0_dp, 4.0_dp, -2.0_dp, 5.0_dp])
   call check('rowMins selected', row_mins(xn, rows = [4, 1, 4], cols = [5, 2, 3], na_rm = .true.), &
      [real(dp) :: 5.0_dp, 3.0_dp, 5.0_dp])
   call check('rowMins missing', row_mins(xn), &
      [real(dp) :: nan, 4.0_dp, nan, 5.0_dp])
   call check('rowMaxs', row_maxs(x), &
      [real(dp) :: 5.0_dp, 9.0_dp, 2.0_dp, 9.0_dp])
   call check('rowMaxs selected', row_maxs(xn, rows = [4, 1, 4], cols = [5, 2, 3], na_rm = .true.), &
      [real(dp) :: 8.0_dp, 5.0_dp, 8.0_dp])
   call check('rowMaxs missing', row_maxs(xn), &
      [real(dp) :: nan, 9.0_dp, nan, 9.0_dp])
   call check('rowLogSumExps', row_log_sum_exps(x), &
      [real(dp) :: 5.4519143959375933_dp, 9.5591784216424749_dp, 2.4519143959375933_dp, 9.4519143959375942_dp])
   call check('rowLogSumExps selected', row_log_sum_exps(xn, rows = [4, 1, 4], cols = [5, 2, 3], na_rm = .true.), &
      [real(dp) :: 8.3490122167681857_dp, 5.1269280110429722_dp, 8.3490122167681857_dp])
   call check('rowLogSumExps missing', row_log_sum_exps(xn), &
      [real(dp) :: nan, 9.5591784216424749_dp, nan, 9.4519143959375942_dp])
   call check('rowCounts', row_counts(x), &
      [real(dp) :: 1.0_dp, 0.0_dp, 1.0_dp, 0.0_dp])
   call check('rowCounts selected', row_counts(xn, rows = [4, 1, 4], cols = [5, 2, 3], na_rm = .true.), &
      [real(dp) :: 0.0_dp, 0.0_dp, 0.0_dp])
   call check('rowCounts missing', row_counts(xn), &
      [real(dp) :: nan, 0.0_dp, nan, 0.0_dp])
   call check('rowAlls', row_alls(x), &
      [real(dp) :: 1.0_dp, 1.0_dp, 0.0_dp, 1.0_dp])
   call check('rowAlls selected', row_alls(xn, rows = [4, 1, 4], cols = [5, 2, 3], na_rm = .true.), &
      [real(dp) :: 1.0_dp, 1.0_dp, 1.0_dp])
   call check('rowAlls missing', row_alls(xn), &
      [real(dp) :: nan, 1.0_dp, 0.0_dp, 1.0_dp])
   call check('rowAnys', row_anys(x), &
      [real(dp) :: 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp])
   call check('rowAnys selected', row_anys(xn, rows = [4, 1, 4], cols = [5, 2, 3], na_rm = .true.), &
      [real(dp) :: 1.0_dp, 1.0_dp, 1.0_dp])
   call check('rowAnys missing', row_anys(xn), &
      [real(dp) :: 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp])
   call check('rowProds', row_prods(x), &
      [real(dp) :: 120.0_dp, 9216.0_dp, 0.0_dp, 15120.0_dp])
   call check('rowProds selected', row_prods(xn, rows = [4, 1, 4], cols = [5, 2, 3], na_rm = .true.), &
      [real(dp) :: 280.0_dp, 15.0_dp, 280.0_dp])
   call check('rowProds missing', row_prods(xn), &
      [real(dp) :: nan, 9216.0_dp, nan, 15120.0_dp])
   call check_matrix('rowRanges', row_ranges(x), &
      [real(dp) :: 1.0_dp, 4.0_dp, -2.0_dp, 5.0_dp, &
      5.0_dp, 9.0_dp, 2.0_dp, 9.0_dp])
   call check_matrix('rowRanges missing', row_ranges(xn), &
      [real(dp) :: nan, 4.0_dp, nan, 5.0_dp, &
      nan, 9.0_dp, nan, 9.0_dp])
   call check_matrix('rowQuantiles', row_quantiles(x, probs = [0.0_dp, 0.25_dp, 0.5_dp, 0.75_dp, 1.0_dp]), &
      [real(dp) :: 1.0_dp, 4.0_dp, -2.0_dp, 5.0_dp, &
      2.0_dp, 4.0_dp, -1.0_dp, 6.0_dp, &
      3.0_dp, 8.0_dp, 0.0_dp, 7.0_dp, &
      4.0_dp, 8.0_dp, 1.0_dp, 8.0_dp, &
      5.0_dp, 9.0_dp, 2.0_dp, 9.0_dp])
   call check_matrix('rowRanks', row_ranks(x, ties_method = 'average'), &
      [real(dp) :: 1.0_dp, 1.5_dp, 1.0_dp, 5.0_dp, &
      2.0_dp, 1.5_dp, 2.0_dp, 4.0_dp, &
      3.0_dp, 3.5_dp, 3.0_dp, 3.0_dp, &
      4.0_dp, 3.5_dp, 4.0_dp, 2.0_dp, &
      5.0_dp, 5.0_dp, 5.0_dp, 1.0_dp])
   call check_matrix('rowRanks missing', row_ranks(xn, ties_method = 'average'), &
      [real(dp) :: 1.0_dp, 1.5_dp, 1.0_dp, 5.0_dp, &
      nan, 1.5_dp, 2.0_dp, 4.0_dp, &
      2.0_dp, 3.5_dp, 3.0_dp, 3.0_dp, &
      3.0_dp, 3.5_dp, nan, 2.0_dp, &
      4.0_dp, 5.0_dp, 4.0_dp, 1.0_dp])
   call check_matrix('rowCumsums', row_cumsums(x), &
      [real(dp) :: 1.0_dp, 4.0_dp, -2.0_dp, 9.0_dp, &
      3.0_dp, 8.0_dp, -3.0_dp, 17.0_dp, &
      6.0_dp, 16.0_dp, -3.0_dp, 24.0_dp, &
      10.0_dp, 24.0_dp, -2.0_dp, 30.0_dp, &
      15.0_dp, 33.0_dp, 0.0_dp, 35.0_dp])
   call check_matrix('rowCumsums missing', row_cumsums(xn), &
      [real(dp) :: 1.0_dp, 4.0_dp, -2.0_dp, 9.0_dp, &
      nan, 8.0_dp, -3.0_dp, 17.0_dp, &
      nan, 16.0_dp, -3.0_dp, 24.0_dp, &
      nan, 24.0_dp, nan, 30.0_dp, &
      nan, 33.0_dp, nan, 35.0_dp])
   call check_matrix('rowCumprods', row_cumprods(x), &
      [real(dp) :: 1.0_dp, 4.0_dp, -2.0_dp, 9.0_dp, &
      2.0_dp, 16.0_dp, 2.0_dp, 72.0_dp, &
      6.0_dp, 128.0_dp, 0.0_dp, 504.0_dp, &
      24.0_dp, 1024.0_dp, 0.0_dp, 3024.0_dp, &
      120.0_dp, 9216.0_dp, 0.0_dp, 15120.0_dp])
   call check_matrix('rowCumprods missing', row_cumprods(xn), &
      [real(dp) :: 1.0_dp, 4.0_dp, -2.0_dp, 9.0_dp, &
      nan, 16.0_dp, 2.0_dp, 72.0_dp, &
      nan, 128.0_dp, 0.0_dp, 504.0_dp, &
      nan, 1024.0_dp, nan, 3024.0_dp, &
      nan, 9216.0_dp, nan, 15120.0_dp])
   call check_matrix('rowCummins', row_cummins(x), &
      [real(dp) :: 1.0_dp, 4.0_dp, -2.0_dp, 9.0_dp, &
      1.0_dp, 4.0_dp, -2.0_dp, 8.0_dp, &
      1.0_dp, 4.0_dp, -2.0_dp, 7.0_dp, &
      1.0_dp, 4.0_dp, -2.0_dp, 6.0_dp, &
      1.0_dp, 4.0_dp, -2.0_dp, 5.0_dp])
   call check_matrix('rowCummins missing', row_cummins(xn), &
      [real(dp) :: 1.0_dp, 4.0_dp, -2.0_dp, 9.0_dp, &
      nan, 4.0_dp, -2.0_dp, 8.0_dp, &
      nan, 4.0_dp, -2.0_dp, 7.0_dp, &
      nan, 4.0_dp, nan, 6.0_dp, &
      nan, 4.0_dp, nan, 5.0_dp])
   call check_matrix('rowCummaxs', row_cummaxs(x), &
      [real(dp) :: 1.0_dp, 4.0_dp, -2.0_dp, 9.0_dp, &
      2.0_dp, 4.0_dp, -1.0_dp, 9.0_dp, &
      3.0_dp, 8.0_dp, 0.0_dp, 9.0_dp, &
      4.0_dp, 8.0_dp, 1.0_dp, 9.0_dp, &
      5.0_dp, 9.0_dp, 2.0_dp, 9.0_dp])
   call check_matrix('rowCummaxs missing', row_cummaxs(xn), &
      [real(dp) :: 1.0_dp, 4.0_dp, -2.0_dp, 9.0_dp, &
      nan, 4.0_dp, -1.0_dp, 9.0_dp, &
      nan, 8.0_dp, 0.0_dp, 9.0_dp, &
      nan, 8.0_dp, nan, 9.0_dp, &
      nan, 9.0_dp, nan, 9.0_dp])
   call check_matrix('rowDiffs', row_diffs(x), &
      [real(dp) :: 1.0_dp, 0.0_dp, 1.0_dp, -1.0_dp, &
      1.0_dp, 4.0_dp, 1.0_dp, -1.0_dp, &
      1.0_dp, 0.0_dp, 1.0_dp, -1.0_dp, &
      1.0_dp, 1.0_dp, 1.0_dp, -1.0_dp])
   call check_matrix('rowDiffs missing', row_diffs(xn), &
      [real(dp) :: nan, 0.0_dp, 1.0_dp, -1.0_dp, &
      nan, 4.0_dp, 1.0_dp, -1.0_dp, &
      1.0_dp, 0.0_dp, nan, -1.0_dp, &
      1.0_dp, 1.0_dp, nan, -1.0_dp])
   call check('rowOrderStats', row_order_stats(x, 2), &
      [real(dp) :: 2.0_dp, 4.0_dp, -1.0_dp, 6.0_dp])
   call check('rowAnyMissings', merge(1.0_dp, 0.0_dp, row_any_missings(xn)), &
      [real(dp) :: 1.0_dp, 0.0_dp, 1.0_dp, 0.0_dp])
   call check('rowAnyNAs', merge(1.0_dp, 0.0_dp, row_any_nas(xn)), &
      [real(dp) :: 1.0_dp, 0.0_dp, 1.0_dp, 0.0_dp])
   call check('rowVarDiffs', row_var_diffs(x, diff = 2), &
      [real(dp) :: 0.0_dp, 4.083333333333333_dp, 0.0_dp, 0.0_dp])
   call check('rowVarDiffs missing', row_var_diffs(xn, na_rm = .true., trim = 0.2_dp), &
      [real(dp) :: 0.16666666666666666_dp, 1.7916666666666667_dp, 0.16666666666666666_dp, 0.0_dp])
   call check('rowSdDiffs', row_sd_diffs(x, diff = 2), &
      [real(dp) :: 0.0_dp, 2.0207259421636898_dp, 0.0_dp, 0.0_dp])
   call check('rowSdDiffs missing', row_sd_diffs(xn, na_rm = .true., trim = 0.2_dp), &
      [real(dp) :: 0.40824829046386296_dp, 1.3385315336840842_dp, 0.40824829046386296_dp, 0.0_dp])
   call check('rowMadDiffs', row_mad_diffs(x, diff = 2), &
      [real(dp) :: 0.0_dp, 2.2238999999999995_dp, 0.0_dp, 0.0_dp])
   call check('rowMadDiffs missing', row_mad_diffs(xn, na_rm = .true., trim = 0.2_dp), &
      [real(dp) :: 0.0_dp, 0.52417825689358766_dp, 0.0_dp, 0.0_dp])
   call check('rowIQRDiffs', row_iqr_diffs(x, diff = 2), &
      [real(dp) :: 0.0_dp, 1.9999999999999996_dp, 0.0_dp, 0.0_dp])
   call check('rowIQRDiffs missing', row_iqr_diffs(xn, na_rm = .true., trim = 0.2_dp), &
      [real(dp) :: 0.35355339059327373_dp, 1.2374368670764582_dp, 0.35355339059327373_dp, 0.0_dp])
   call check('rowWeightedMeans no weights', row_weighted_means(x), &
      [real(dp) :: 3.0_dp, 6.5999999999999996_dp, 0.0_dp, 7.0_dp])
   call check('rowWeightedMeans', row_weighted_means(x, wr), &
      [real(dp) :: 3.4285714285714284_dp, 7.5714285714285712_dp, 0.42857142857142855_dp, 6.5714285714285712_dp])
   call check('rowWeightedMeans missing', row_weighted_means(xn, wr, rows = [4, 2], cols = [5, 2, 1], na_rm = .true.), &
      [real(dp) :: 7.0_dp, 6.5_dp])
   call check('rowWeightedMedians no weights', row_weighted_medians(x), &
      [real(dp) :: 3.0_dp, 8.0_dp, 0.0_dp, 7.0_dp])
   call check('rowWeightedMedians', row_weighted_medians(x, wr), &
      [real(dp) :: 3.6000000000000001_dp, 8.0_dp, 0.60000000000000009_dp, 6.4000000000000004_dp])
   call check('rowWeightedMedians missing', row_weighted_medians(xn, wr, rows = [4, 2], cols = [5, 2, 1], na_rm = .true.), &
      [real(dp) :: 7.0_dp, 6.5_dp])
   call check('rowWeightedVars no weights', row_weighted_vars(x), &
      [real(dp) :: 2.5_dp, 5.7999999999999998_dp, 2.5_dp, 2.5_dp])
   call check('rowWeightedVars', row_weighted_vars(x, wr), &
      [real(dp) :: 1.6190476190476191_dp, 2.6190476190476191_dp, 1.6190476190476188_dp, 1.6190476190476188_dp])
   call check('rowWeightedVars missing', row_weighted_vars(xn, wr, rows = [4, 2], cols = [5, 2, 1], na_rm = .true.), &
      [real(dp) :: 8.0_dp, 12.5_dp])
   call check('rowWeightedSds no weights', row_weighted_sds(x), &
      [real(dp) :: 1.5811388300841898_dp, 2.4083189157584592_dp, 1.5811388300841898_dp, 1.5811388300841898_dp])
   call check('rowWeightedSds', row_weighted_sds(x, wr), &
      [real(dp) :: 1.2724180205607036_dp, 1.6183471874253741_dp, 1.2724180205607034_dp, 1.2724180205607034_dp])
   call check('rowWeightedSds missing', row_weighted_sds(xn, wr, rows = [4, 2], cols = [5, 2, 1], na_rm = .true.), &
      [real(dp) :: 2.8284271247461903_dp, 3.5355339059327378_dp])
   call check('rowWeightedMads no weights', row_weighted_mads(x), &
      [real(dp) :: 1.4825999999999999_dp, 1.4825999999999999_dp, 1.4825999999999999_dp, 1.4825999999999999_dp])
   call check('rowWeightedMads', row_weighted_mads(x, wr), &
      [real(dp) :: 0.83025599999999999_dp, 6.584066625237028e-16_dp, 0.83025599999999999_dp, 0.83025599999999966_dp])
   call check('rowWeightedMads missing', row_weighted_mads(xn, wr, rows = [4, 2], cols = [5, 2, 1], na_rm = .true.), &
      [real(dp) :: 2.9651999999999998_dp, 3.7064999999999997_dp])
   call check_matrix('rowTabulates', real(row_tabulates(xn, [1.0_dp, 4.0_dp, nan]), dp), &
      [real(dp) :: 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
      1.0_dp, 2.0_dp, 0.0_dp, 0.0_dp, &
      1.0_dp, 0.0_dp, 1.0_dp, 0.0_dp])
   call check('rowCollapse', row_collapse(x, [3, 1, 2]), &
      [real(dp) :: 3.0_dp, 4.0_dp, -1.0_dp, 7.0_dp])
   call check('colSums2', col_sums2(x), &
      [real(dp) :: 12.0_dp, 13.0_dp, 18.0_dp, 19.0_dp, &
      21.0_dp])
   call check('colSums2 selected', col_sums2(xn, rows = [4, 1, 4], cols = [5, 2, 3], na_rm = .true.), &
      [real(dp) :: 15.0_dp, 16.0_dp, 17.0_dp])
   call check('colSums2 missing', col_sums2(xn), &
      [real(dp) :: 12.0_dp, nan, 18.0_dp, nan, &
      21.0_dp])
   call check('colMeans2', col_means2(x), &
      [real(dp) :: 3.0_dp, 3.25_dp, 4.5_dp, 4.75_dp, &
      5.25_dp])
   call check('colMeans2 selected', col_means2(xn, rows = [4, 1, 4], cols = [5, 2, 3], na_rm = .true.), &
      [real(dp) :: 5.0_dp, 8.0_dp, 5.666666666666667_dp])
   call check('colMeans2 missing', col_means2(xn), &
      [real(dp) :: 3.0_dp, nan, 4.5_dp, nan, &
      5.25_dp])
   call check('colVars', col_vars(x), &
      [real(dp) :: 22.0_dp, 14.25_dp, 13.666666666666666_dp, 8.9166666666666661_dp, &
      8.25_dp])
   call check('colVars selected', col_vars(xn, rows = [4, 1, 4], cols = [5, 2, 3], na_rm = .true.), &
      [real(dp) :: 0.0_dp, 0.0_dp, 5.333333333333333_dp])
   call check('colVars missing', col_vars(xn), &
      [real(dp) :: 22.0_dp, nan, 13.666666666666666_dp, nan, &
      8.25_dp])
   call check('colSds', col_sds(x), &
      [real(dp) :: 4.6904157598234297_dp, 3.7749172176353749_dp, 3.6968455021364721_dp, 2.9860788111948193_dp, &
      2.8722813232690143_dp])
   call check('colSds selected', col_sds(xn, rows = [4, 1, 4], cols = [5, 2, 3], na_rm = .true.), &
      [real(dp) :: 0.0_dp, 0.0_dp, 2.3094010767585029_dp])
   call check('colSds missing', col_sds(xn), &
      [real(dp) :: 4.6904157598234297_dp, nan, 3.6968455021364721_dp, nan, &
      2.8722813232690143_dp])
   call check('colMedians', col_medians(x), &
      [real(dp) :: 2.5_dp, 3.0_dp, 5.0_dp, 5.0_dp, &
      5.0_dp])
   call check('colMedians selected', col_medians(xn, rows = [4, 1, 4], cols = [5, 2, 3], na_rm = .true.), &
      [real(dp) :: 5.0_dp, 8.0_dp, 7.0_dp])
   call check('colMedians missing', col_medians(xn), &
      [real(dp) :: 2.5_dp, nan, 5.0_dp, nan, &
      5.0_dp])
   call check('colMads', col_mads(x), &
      [real(dp) :: 4.4478_dp, 3.7064999999999997_dp, 3.7064999999999997_dp, 2.9651999999999998_dp, &
      2.2239_dp])
   call check('colMads selected', col_mads(xn, rows = [4, 1, 4], cols = [5, 2, 3], na_rm = .true.), &
      [real(dp) :: 0.0_dp, 0.0_dp, 0.0_dp])
   call check('colMads missing', col_mads(xn), &
      [real(dp) :: 4.4478_dp, nan, 3.7064999999999997_dp, nan, &
      2.2239_dp])
   call check('colIQRs', col_iqrs(x), &
      [real(dp) :: 5.0_dp, 3.75_dp, 5.0_dp, 3.25_dp, &
      1.75_dp])
   call check('colIQRs selected', col_iqrs(xn, rows = [4, 1, 4], cols = [5, 2, 3], na_rm = .true.), &
      [real(dp) :: 0.0_dp, 0.0_dp, 2.0_dp])
   call check('colIQRs missing', col_iqrs(xn), &
      [real(dp) :: 5.0_dp, nan, 5.0_dp, nan, &
      1.75_dp])
   call check('colMins', col_mins(x), &
      [real(dp) :: -2.0_dp, -1.0_dp, 0.0_dp, 1.0_dp, &
      2.0_dp])
   call check('colMins selected', col_mins(xn, rows = [4, 1, 4], cols = [5, 2, 3], na_rm = .true.), &
      [real(dp) :: 5.0_dp, 8.0_dp, 3.0_dp])
   call check('colMins missing', col_mins(xn), &
      [real(dp) :: -2.0_dp, nan, 0.0_dp, nan, &
      2.0_dp])
   call check('colMaxs', col_maxs(x), &
      [real(dp) :: 9.0_dp, 8.0_dp, 8.0_dp, 8.0_dp, &
      9.0_dp])
   call check('colMaxs selected', col_maxs(xn, rows = [4, 1, 4], cols = [5, 2, 3], na_rm = .true.), &
      [real(dp) :: 5.0_dp, 8.0_dp, 7.0_dp])
   call check('colMaxs missing', col_maxs(xn), &
      [real(dp) :: 9.0_dp, nan, 8.0_dp, nan, &
      9.0_dp])
   call check('colLogSumExps', col_log_sum_exps(x), &
      [real(dp) :: 9.0070650946661335_dp, 8.0207020274885501_dp, 8.3184194402021934_dp, 8.1437217477186152_dp, &
      9.036855572009781_dp])
   call check('colLogSumExps selected', col_log_sum_exps(xn, rows = [4, 1, 4], cols = [5, 2, 3], na_rm = .true.), &
      [real(dp) :: 6.0986122886681091_dp, 8.6931471805599454_dp, 7.7022633214390952_dp])
   call check('colLogSumExps missing', col_log_sum_exps(xn), &
      [real(dp) :: 9.0070650946661335_dp, nan, 8.3184194402021934_dp, nan, &
      9.036855572009781_dp])
   call check('colCounts', col_counts(x), &
      [real(dp) :: 1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, &
      0.0_dp])
   call check('colCounts selected', col_counts(xn, rows = [4, 1, 4], cols = [5, 2, 3], na_rm = .true.), &
      [real(dp) :: 0.0_dp, 0.0_dp, 0.0_dp])
   call check('colCounts missing', col_counts(xn), &
      [real(dp) :: 1.0_dp, nan, 0.0_dp, nan, &
      0.0_dp])
   call check('colAlls', col_alls(x), &
      [real(dp) :: 1.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, &
      1.0_dp])
   call check('colAlls selected', col_alls(xn, rows = [4, 1, 4], cols = [5, 2, 3], na_rm = .true.), &
      [real(dp) :: 1.0_dp, 1.0_dp, 1.0_dp])
   call check('colAlls missing', col_alls(xn), &
      [real(dp) :: 1.0_dp, nan, 0.0_dp, nan, &
      1.0_dp])
   call check('colAnys', col_anys(x), &
      [real(dp) :: 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, &
      1.0_dp])
   call check('colAnys selected', col_anys(xn, rows = [4, 1, 4], cols = [5, 2, 3], na_rm = .true.), &
      [real(dp) :: 1.0_dp, 1.0_dp, 1.0_dp])
   call check('colAnys missing', col_anys(xn), &
      [real(dp) :: 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, &
      1.0_dp])
   call check('colProds', col_prods(x), &
      [real(dp) :: -72.0_dp, -64.0_dp, 0.0_dp, 192.0_dp, &
      450.0_dp])
   call check('colProds selected', col_prods(xn, rows = [4, 1, 4], cols = [5, 2, 3], na_rm = .true.), &
      [real(dp) :: 125.0_dp, 64.0_dp, 147.0_dp])
   call check('colProds missing', col_prods(xn), &
      [real(dp) :: -72.0_dp, nan, 0.0_dp, nan, &
      450.0_dp])
   call check_matrix('colRanges', col_ranges(x), &
      [real(dp) :: -2.0_dp, -1.0_dp, 0.0_dp, 1.0_dp, &
      2.0_dp, 9.0_dp, 8.0_dp, 8.0_dp, &
      8.0_dp, 9.0_dp])
   call check_matrix('colRanges missing', col_ranges(xn), &
      [real(dp) :: -2.0_dp, nan, 0.0_dp, nan, &
      2.0_dp, 9.0_dp, nan, 8.0_dp, &
      nan, 9.0_dp])
   call check_matrix('colQuantiles', col_quantiles(x, probs = [0.0_dp, 0.25_dp, 0.5_dp, 0.75_dp, 1.0_dp]), &
      [real(dp) :: -2.0_dp, -1.0_dp, 0.0_dp, 1.0_dp, &
      2.0_dp, 0.25_dp, 1.25_dp, 2.25_dp, &
      3.25_dp, 4.25_dp, 2.5_dp, 3.0_dp, &
      5.0_dp, 5.0_dp, 5.0_dp, 5.25_dp, &
      5.0_dp, 7.25_dp, 6.5_dp, 6.0_dp, &
      9.0_dp, 8.0_dp, 8.0_dp, 8.0_dp, &
      9.0_dp])
   call check_matrix('colRanks', col_ranks(x, ties_method = 'average'), &
      [real(dp) :: 2.0_dp, 2.0_dp, 2.0_dp, 2.0_dp, &
      2.5_dp, 3.0_dp, 3.0_dp, 4.0_dp, &
      4.0_dp, 4.0_dp, 1.0_dp, 1.0_dp, &
      1.0_dp, 1.0_dp, 1.0_dp, 4.0_dp, &
      4.0_dp, 3.0_dp, 3.0_dp, 2.5_dp])
   call check_matrix('colRanks missing', col_ranks(xn, ties_method = 'average'), &
      [real(dp) :: 2.0_dp, nan, 2.0_dp, 1.0_dp, &
      2.5_dp, 3.0_dp, 2.0_dp, 4.0_dp, &
      3.0_dp, 4.0_dp, 1.0_dp, 1.0_dp, &
      1.0_dp, nan, 1.0_dp, 4.0_dp, &
      3.0_dp, 3.0_dp, 2.0_dp, 2.5_dp])
   call check_matrix('colCumsums', col_cumsums(x), &
      [real(dp) :: 1.0_dp, 5.0_dp, 3.0_dp, 12.0_dp, &
      2.0_dp, 6.0_dp, 5.0_dp, 13.0_dp, &
      3.0_dp, 11.0_dp, 11.0_dp, 18.0_dp, &
      4.0_dp, 12.0_dp, 13.0_dp, 19.0_dp, &
      5.0_dp, 14.0_dp, 16.0_dp, 21.0_dp])
   call check_matrix('colCumsums missing', col_cumsums(xn), &
      [real(dp) :: 1.0_dp, 5.0_dp, 3.0_dp, 12.0_dp, &
      nan, nan, nan, nan, &
      3.0_dp, 11.0_dp, 11.0_dp, 18.0_dp, &
      4.0_dp, 12.0_dp, nan, nan, &
      5.0_dp, 14.0_dp, 16.0_dp, 21.0_dp])
   call check_matrix('colCumprods', col_cumprods(x), &
      [real(dp) :: 1.0_dp, 4.0_dp, -8.0_dp, -72.0_dp, &
      2.0_dp, 8.0_dp, -8.0_dp, -64.0_dp, &
      3.0_dp, 24.0_dp, 0.0_dp, 0.0_dp, &
      4.0_dp, 32.0_dp, 32.0_dp, 192.0_dp, &
      5.0_dp, 45.0_dp, 90.0_dp, 450.0_dp])
   call check_matrix('colCumprods missing', col_cumprods(xn), &
      [real(dp) :: 1.0_dp, 4.0_dp, -8.0_dp, -72.0_dp, &
      nan, nan, nan, nan, &
      3.0_dp, 24.0_dp, 0.0_dp, 0.0_dp, &
      4.0_dp, 32.0_dp, nan, nan, &
      5.0_dp, 45.0_dp, 90.0_dp, 450.0_dp])
   call check_matrix('colCummins', col_cummins(x), &
      [real(dp) :: 1.0_dp, 1.0_dp, -2.0_dp, -2.0_dp, &
      2.0_dp, 2.0_dp, -1.0_dp, -1.0_dp, &
      3.0_dp, 3.0_dp, 0.0_dp, 0.0_dp, &
      4.0_dp, 4.0_dp, 1.0_dp, 1.0_dp, &
      5.0_dp, 5.0_dp, 2.0_dp, 2.0_dp])
   call check_matrix('colCummins missing', col_cummins(xn), &
      [real(dp) :: 1.0_dp, 1.0_dp, -2.0_dp, -2.0_dp, &
      nan, nan, nan, nan, &
      3.0_dp, 3.0_dp, 0.0_dp, 0.0_dp, &
      4.0_dp, 4.0_dp, nan, nan, &
      5.0_dp, 5.0_dp, 2.0_dp, 2.0_dp])
   call check_matrix('colCummaxs', col_cummaxs(x), &
      [real(dp) :: 1.0_dp, 4.0_dp, 4.0_dp, 9.0_dp, &
      2.0_dp, 4.0_dp, 4.0_dp, 8.0_dp, &
      3.0_dp, 8.0_dp, 8.0_dp, 8.0_dp, &
      4.0_dp, 8.0_dp, 8.0_dp, 8.0_dp, &
      5.0_dp, 9.0_dp, 9.0_dp, 9.0_dp])
   call check_matrix('colCummaxs missing', col_cummaxs(xn), &
      [real(dp) :: 1.0_dp, 4.0_dp, 4.0_dp, 9.0_dp, &
      nan, nan, nan, nan, &
      3.0_dp, 8.0_dp, 8.0_dp, 8.0_dp, &
      4.0_dp, 8.0_dp, nan, nan, &
      5.0_dp, 9.0_dp, 9.0_dp, 9.0_dp])
   call check_matrix('colDiffs', col_diffs(x), &
      [real(dp) :: 3.0_dp, -6.0_dp, 11.0_dp, 2.0_dp, &
      -5.0_dp, 9.0_dp, 5.0_dp, -8.0_dp, &
      7.0_dp, 4.0_dp, -7.0_dp, 5.0_dp, &
      4.0_dp, -7.0_dp, 3.0_dp])
   call check_matrix('colDiffs missing', col_diffs(xn), &
      [real(dp) :: 3.0_dp, -6.0_dp, 11.0_dp, nan, &
      -5.0_dp, 9.0_dp, 5.0_dp, -8.0_dp, &
      7.0_dp, 4.0_dp, nan, nan, &
      4.0_dp, -7.0_dp, 3.0_dp])
   call check('colOrderStats', col_order_stats(x, 2), &
      [real(dp) :: 1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, &
      5.0_dp])
   call check('colAnyMissings', merge(1.0_dp, 0.0_dp, col_any_missings(xn)), &
      [real(dp) :: 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, &
      0.0_dp])
   call check('colAnyNAs', merge(1.0_dp, 0.0_dp, col_any_nas(xn)), &
      [real(dp) :: 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, &
      0.0_dp])
   call check('colVarDiffs', col_var_diffs(x, diff = 2), &
      [real(dp) :: 84.5_dp, 55.125_dp, 98.0_dp, 66.125_dp, &
      55.125_dp])
   call check('colVarDiffs missing', col_var_diffs(xn, na_rm = .true., trim = 0.2_dp), &
      [real(dp) :: 36.166666666666664_dp, 49.0_dp, 33.166666666666664_dp, 9.0_dp, &
      18.5_dp])
   call check('colSdDiffs', col_sd_diffs(x, diff = 2), &
      [real(dp) :: 9.1923881554251157_dp, 7.4246212024587468_dp, 9.8994949366116636_dp, 8.1317279836452947_dp, &
      7.4246212024587468_dp])
   call check('colSdDiffs missing', col_sd_diffs(xn, na_rm = .true., trim = 0.2_dp), &
      [real(dp) :: 6.0138728508895705_dp, 6.9999999999999991_dp, 5.7590508477236648_dp, 2.9999999999999996_dp, &
      4.3011626335213125_dp])
   call check('colMadDiffs', col_mad_diffs(x, diff = 2), &
      [real(dp) :: 9.6368999999999971_dp, 7.783649999999998_dp, 10.378199999999998_dp, 8.5249499999999969_dp, &
      7.783649999999998_dp])
   call check('colMadDiffs missing', col_mad_diffs(xn, na_rm = .true., trim = 0.2_dp), &
      [real(dp) :: 8.3868521102974025_dp, 7.338495596510227_dp, 2.0967130275743506_dp, 3.1450695413615257_dp, &
      1.0483565137871753_dp])
   call check('colIQRDiffs', col_iqr_diffs(x, diff = 2), &
      [real(dp) :: 6.4999999999999982_dp, 5.2499999999999991_dp, 6.9999999999999982_dp, 5.7499999999999991_dp, &
      5.2499999999999991_dp])
   call check('colIQRDiffs missing', col_iqr_diffs(xn, na_rm = .true., trim = 0.2_dp), &
      [real(dp) :: 6.0104076400856536_dp, 4.9497474683058327_dp, 5.3033008588991057_dp, 2.1213203435596424_dp, &
      3.8890872965260113_dp])
   call check('colWeightedMeans no weights', col_weighted_means(x), &
      [real(dp) :: 3.0_dp, 3.25_dp, 4.5_dp, 4.75_dp, &
      5.25_dp])
   call check('colWeightedMeans', col_weighted_means(x, wc), &
      [real(dp) :: 6.0_dp, 5.666666666666667_dp, 6.6666666666666661_dp, 6.333333333333333_dp, &
      6.333333333333333_dp])
   call check('colWeightedMeans missing', col_weighted_means(xn, wc, rows = [4, 2], cols = [5, 2, 1], na_rm = .true.), &
      [real(dp) :: 6.5999999999999996_dp, 6.4000000000000004_dp, 7.0_dp])
   call check('colWeightedMedians no weights', col_weighted_medians(x), &
      [real(dp) :: 2.5_dp, 3.0_dp, 5.0_dp, 5.0_dp, &
      5.0_dp])
   call check('colWeightedMedians', col_weighted_medians(x, wc), &
      [real(dp) :: 6.0_dp, 5.5999999999999996_dp, 7.2000000000000002_dp, 6.4000000000000004_dp, &
      5.8000000000000007_dp])
   call check('colWeightedMedians missing', col_weighted_medians(xn, wc, rows = [4, 2], cols = [5, 2, 1], na_rm = .true.), &
      [real(dp) :: 6.5999999999999996_dp, 6.4000000000000004_dp, 7.0_dp])
   call check('colWeightedVars no weights', col_weighted_vars(x), &
      [real(dp) :: 22.0_dp, 14.25_dp, 13.666666666666666_dp, 8.9166666666666661_dp, &
      8.25_dp])
   call check('colWeightedVars', col_weighted_vars(x, wc), &
      [real(dp) :: 12.0_dp, 7.0666666666666664_dp, 3.4666666666666668_dp, 2.2666666666666666_dp, &
      4.2666666666666666_dp])
   call check('colWeightedVars missing', col_weighted_vars(xn, wc, rows = [4, 2], cols = [5, 2, 1], na_rm = .true.), &
      [real(dp) :: 4.7999999999999998_dp, 4.7999999999999998_dp, 7.5_dp])
   call check('colWeightedSds no weights', col_weighted_sds(x), &
      [real(dp) :: 4.6904157598234297_dp, 3.7749172176353749_dp, 3.6968455021364721_dp, 2.9860788111948193_dp, &
      2.8722813232690143_dp])
   call check('colWeightedSds', col_weighted_sds(x, wc), &
      [real(dp) :: 3.4641016151377544_dp, 2.6583202716502514_dp, 1.8618986725025255_dp, 1.505545305418162_dp, &
      2.0655911179772888_dp])
   call check('colWeightedSds missing', col_weighted_sds(xn, wc, rows = [4, 2], cols = [5, 2, 1], na_rm = .true.), &
      [real(dp) :: 2.1908902300206643_dp, 2.1908902300206643_dp, 2.7386127875258306_dp])
   call check('colWeightedMads no weights', col_weighted_mads(x), &
      [real(dp) :: 4.4478_dp, 3.7064999999999997_dp, 3.7064999999999997_dp, 2.9651999999999998_dp, &
      2.2239_dp])
   call check('colWeightedMads', col_weighted_mads(x, wc), &
      [real(dp) :: 4.1512799999999999_dp, 3.321024_dp, 0.83025599999999999_dp, 1.660512_dp, &
      1.897728000000001_dp])
   call check('colWeightedMads missing', col_weighted_mads(xn, wc, rows = [4, 2], cols = [5, 2, 1], na_rm = .true.), &
      [real(dp) :: 2.8465919999999998_dp, 2.8465919999999998_dp, 3.5582399999999996_dp])
   call check_matrix('colTabulates', real(col_tabulates(xn, [1.0_dp, 4.0_dp, nan]), dp), &
      [real(dp) :: 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
      0.0_dp, 1.0_dp, 1.0_dp, 0.0_dp, &
      1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, &
      0.0_dp, 1.0_dp, 0.0_dp])
   call check('colCollapse', col_collapse(x, [3, 1, 2]), &
      [real(dp) :: -2.0_dp, 2.0_dp, 8.0_dp, 1.0_dp, &
      5.0_dp])
   call check('sum2', [sum2(v)], &
      [real(dp) :: 24.0_dp])
   call check('sum2 subset', [sum2(vn, idxs = [6, 1, 4, 2], na_rm = .true.)], &
      [real(dp) :: 7.0_dp])
   call check('mean2', [mean2(v)], &
      [real(dp) :: 4.0_dp])
   call check('mean2 subset', [mean2(vn, idxs = [6, 1, 4, 2], na_rm = .true.)], &
      [real(dp) :: 2.3333333333333335_dp])
   call check('product', [product(v)], &
      [real(dp) :: 504.00000000000006_dp])
   call check('product subset', [product(vn, idxs = [6, 1, 4, 2], na_rm = .true.)], &
      [real(dp) :: 7.9999999999999982_dp])
   call check('allValue', [all_value(v)], &
      [real(dp) :: 1.0_dp])
   call check('allValue subset', [all_value(vn, idxs = [6, 1, 4, 2], na_rm = .true.)], &
      [real(dp) :: 1.0_dp])
   call check('anyValue', [any_value(v)], &
      [real(dp) :: 1.0_dp])
   call check('anyValue subset', [any_value(vn, idxs = [6, 1, 4, 2], na_rm = .true.)], &
      [real(dp) :: 1.0_dp])
   call check('count', [count_value(v)], &
      [real(dp) :: 2.0_dp])
   call check('count subset', [count_value(vn, idxs = [6, 1, 4, 2], na_rm = .true.)], &
      [real(dp) :: 1.0_dp])
   call check('logSumExp', [log_sum_exp(v)], &
      [real(dp) :: 9.1342301811898388_dp])
   call check('logSumExp subset', [log_sum_exp(vn, idxs = [6, 1, 4, 2], na_rm = .true.)], &
      [real(dp) :: 4.1698460195562861_dp])
   call check('iqr', [iqr(v)], &
      [real(dp) :: 5.0_dp])
   call check('iqr subset', [iqr(vn, idxs = [6, 1, 4, 2], na_rm = .true.)], &
      [real(dp) :: 1.5_dp])
   call check('anyMissing', [merge(1.0_dp, 0.0_dp, any_missing(vn))], &
      [real(dp) :: 1.0_dp])
   call check('diff2', diff2(v, lag = 2, differences = 2), &
      [real(dp) :: -1.0_dp, 1.0_dp])
   call check('varDiff', [var_diff(v, diff = 2, trim = 0.2_dp)], &
      [real(dp) :: 53.5_dp])
   call check('sdDiff', [sd_diff(v, diff = 2, trim = 0.2_dp)], &
      [real(dp) :: 7.3143694191638948_dp])
   call check('madDiff', [mad_diff(v, diff = 2, trim = 0.2_dp)], &
      [real(dp) :: 8.8955999999999982_dp])
   call check('iqrDiff', [iqr_diff(v, diff = 2, trim = 0.2_dp)], &
      [real(dp) :: 11.499999999999998_dp])
   call check('weightedMean', [weighted_mean(v, w)], &
      [real(dp) :: 4.083333333333333_dp])
   call check('weightedMean missing subset', [weighted_mean(vn, w, idxs = [6, 3, 2, 1], na_rm = .true.)], &
      [real(dp) :: 2.6666666666666665_dp])
   call check('weightedMedian', [weighted_median(v, w)], &
      [real(dp) :: 2.0000000000000009_dp])
   call check('weightedMedian missing subset', [weighted_median(vn, w, idxs = [6, 3, 2, 1], na_rm = .true.)], &
      [real(dp) :: 2.6666666666666665_dp])
   call check('weightedVar', [weighted_var(v, w)], &
      [real(dp) :: 13.90151515151515_dp])
   call check('weightedVar missing subset', [weighted_var(vn, w, idxs = [6, 3, 2, 1], na_rm = .true.)], &
      [real(dp) :: 1.3333333333333333_dp])
   call check('weightedSd', [weighted_sd(v, w)], &
      [real(dp) :: 3.7284735685686643_dp])
   call check('weightedSd missing subset', [weighted_sd(vn, w, idxs = [6, 3, 2, 1], na_rm = .true.)], &
      [real(dp) :: 1.1547005383792515_dp])
   call check('weightedMad', [weighted_mad(v, w)], &
      [real(dp) :: 1.4826000000000019_dp])
   call check('weightedMad missing subset', [weighted_mad(vn, w, idxs = [6, 3, 2, 1], na_rm = .true.)], &
      [real(dp) :: 1.3178666666666665_dp])
   call check('binCounts', real(bin_counts(v, [0.0_dp, 2.0_dp, 4.0_dp, 8.0_dp, 10.0_dp]), dp), &
      [real(dp) :: 2.0_dp, 1.0_dp, 2.0_dp, 1.0_dp])
   call check('binCounts right', real(bin_counts(v, [0.0_dp, 2.0_dp, 4.0_dp, 8.0_dp, 10.0_dp], right = .true.), dp), &
      [real(dp) :: 3.0_dp, 1.0_dp, 1.0_dp, 1.0_dp])
   call check('binMeans', bm, &
      [real(dp) :: 1.0_dp, 2.0_dp, 5.5_dp, 9.0_dp])
   call check('signTabulate', real(sign_tabulate([ninf, -2.0_dp, 0.0_dp, nan, nan, inf, 1.0_dp]), dp), &
      [real(dp) :: 2.0_dp, 1.0_dp, 2.0_dp, 2.0_dp, &
      1.0_dp, 1.0_dp])
   call check('indexByRow', real(index_by_row([4, 5], [1, 2, 5, 6, 20]), dp), &
      [real(dp) :: 1.0_dp, 5.0_dp, 17.0_dp, 2.0_dp, &
      20.0_dp])
   call check_matrix('rowAvgsPerColSet', row_avgs_per_col_set(x, reshape([1, 3, 5, 2, 4, 0], [3, 2])), &
      [real(dp) :: 3.0_dp, 7.0_dp, 0.0_dp, 7.0_dp, &
      3.0_dp, 6.0_dp, 0.0_dp, 7.0_dp])
   call check_matrix('colAvgsPerRowSet', col_avgs_per_row_set(x, reshape([1, 3, 4, 2, 4, 0], [3, 2])), &
      [real(dp) :: 2.6666666666666665_dp, 6.5_dp, 3.0_dp, 6.0_dp, &
      3.3333333333333335_dp, 7.5_dp, 3.6666666666666665_dp, 7.0_dp, &
      4.0_dp, 7.0_dp])
   call check_matrix('x_OP_y+', x_op_y(xn, [2.0_dp, nan, 4.0_dp], '+', commute = .true., na_rm = .true.), &
      [real(dp) :: 3.0_dp, 4.0_dp, 2.0_dp, 11.0_dp, &
      nan, 8.0_dp, 1.0_dp, 8.0_dp, &
      7.0_dp, 10.0_dp, 0.0_dp, 11.0_dp, &
      6.0_dp, 8.0_dp, 4.0_dp, 8.0_dp, &
      5.0_dp, 13.0_dp, 4.0_dp, 5.0_dp])
   call check_matrix('x_OP_y-', x_op_y(xn, [2.0_dp, nan, 4.0_dp], '-', commute = .true., na_rm = .true.), &
      [real(dp) :: 1.0_dp, nan, 6.0_dp, -7.0_dp, &
      nan, 0.0_dp, 3.0_dp, nan, &
      1.0_dp, -6.0_dp, nan, -3.0_dp, &
      -2.0_dp, nan, nan, -4.0_dp, &
      nan, -5.0_dp, 0.0_dp, nan])
   call check_matrix('x_OP_y*', x_op_y(xn, [2.0_dp, nan, 4.0_dp], '*', commute = .true., na_rm = .true.), &
      [real(dp) :: 2.0_dp, 4.0_dp, -8.0_dp, 18.0_dp, &
      nan, 16.0_dp, -2.0_dp, 8.0_dp, &
      12.0_dp, 16.0_dp, 0.0_dp, 28.0_dp, &
      8.0_dp, 8.0_dp, 4.0_dp, 12.0_dp, &
      5.0_dp, 36.0_dp, 4.0_dp, 5.0_dp])
   call check_matrix('x_OP_y/', x_op_y(xn, [2.0_dp, nan, 4.0_dp], '/', commute = .true., na_rm = .true.), &
      [real(dp) :: 2.0_dp, nan, -2.0_dp, 0.22222222222222221_dp, &
      nan, 1.0_dp, -2.0_dp, nan, &
      1.3333333333333333_dp, 0.25_dp, nan, 0.5714285714285714_dp, &
      0.5_dp, nan, nan, 0.33333333333333331_dp, &
      nan, 0.44444444444444442_dp, 1.0_dp, nan])
   call check_matrix('t_tx_OP_y+', t_tx_op_y(xn, [2.0_dp, nan, 4.0_dp], '+', commute = .true., na_rm = .true.), &
      [real(dp) :: 3.0_dp, 8.0_dp, -2.0_dp, 11.0_dp, &
      nan, 6.0_dp, 3.0_dp, 8.0_dp, &
      7.0_dp, 8.0_dp, 2.0_dp, 11.0_dp, &
      6.0_dp, 12.0_dp, nan, 8.0_dp, &
      5.0_dp, 11.0_dp, 6.0_dp, 5.0_dp])
   call check_matrix('t_tx_OP_y-', t_tx_op_y(xn, [2.0_dp, nan, 4.0_dp], '-', commute = .true., na_rm = .true.), &
      [real(dp) :: 1.0_dp, 0.0_dp, nan, -7.0_dp, &
      nan, -2.0_dp, 5.0_dp, nan, &
      1.0_dp, nan, 2.0_dp, -3.0_dp, &
      -2.0_dp, -4.0_dp, nan, -4.0_dp, &
      nan, -7.0_dp, 2.0_dp, nan])
   call check_matrix('t_tx_OP_y*', t_tx_op_y(xn, [2.0_dp, nan, 4.0_dp], '*', commute = .true., na_rm = .true.), &
      [real(dp) :: 2.0_dp, 16.0_dp, -2.0_dp, 18.0_dp, &
      nan, 8.0_dp, -4.0_dp, 8.0_dp, &
      12.0_dp, 8.0_dp, 0.0_dp, 28.0_dp, &
      8.0_dp, 32.0_dp, nan, 12.0_dp, &
      5.0_dp, 18.0_dp, 8.0_dp, 5.0_dp])
   call check_matrix('t_tx_OP_y/', t_tx_op_y(xn, [2.0_dp, nan, 4.0_dp], '/', commute = .true., na_rm = .true.), &
      [real(dp) :: 2.0_dp, 1.0_dp, nan, 0.22222222222222221_dp, &
      nan, 0.5_dp, -4.0_dp, nan, &
      1.3333333333333333_dp, nan, inf, 0.5714285714285714_dp, &
      0.5_dp, 0.5_dp, nan, 0.33333333333333331_dp, &
      nan, 0.22222222222222221_dp, 2.0_dp, nan])
   call check('rowMeans centered', row_vars(x, center = row_means2(x)), &
      [real(dp) :: 2.5_dp, 5.7999999999999998_dp, 2.5_dp, 2.5_dp])
   call check('rowMAD centered', row_mads(x, center = row_medians(x), constant = 2.0_dp), &
      [real(dp) :: 2.0_dp, 2.0_dp, 2.0_dp, 2.0_dp])
   call check('rownegative subset', row_sums2(x, rows = [-1, -3], cols = [0, 2, 5]), &
      [real(dp) :: 13.0_dp, 13.0_dp])
   call check('colMeans centered', col_vars(x, center = col_means2(x)), &
      [real(dp) :: 22.0_dp, 14.25_dp, 13.666666666666666_dp, 8.9166666666666661_dp, &
      8.25_dp])
   call check('colMAD centered', col_mads(x, center = col_medians(x), constant = 2.0_dp), &
      [real(dp) :: 6.0_dp, 5.0_dp, 5.0_dp, 4.0_dp, &
      3.0_dp])
   call check('colnegative subset', col_sums2(x, rows = [-1, -3], cols = [0, 2, 5]), &
      [real(dp) :: 12.0_dp, 14.0_dp])
   call check('sum2 empty', [sum2([real(dp) ::])], &
      [real(dp) :: 0.0_dp])
   call check('sum2 infinities', [sum2([ninf, inf])], &
      [real(dp) :: nan])
   call check('mean2 empty', [mean2([real(dp) ::])], &
      [real(dp) :: nan])
   call check('mean2 infinities', [mean2([ninf, inf])], &
      [real(dp) :: nan])
   call check('product empty', [product([real(dp) ::])], &
      [real(dp) :: 1.0_dp])
   call check('product infinities', [product([ninf, inf])], &
      [real(dp) :: ninf])
   call check('logSumExp empty', [log_sum_exp([real(dp) ::])], &
      [real(dp) :: ninf])
   call check('logSumExp infinities', [log_sum_exp([ninf, inf])], &
      [real(dp) :: inf])
   call check('iqr empty', [iqr([real(dp) ::])], &
      [real(dp) :: nan])
   call check('iqr infinities', [iqr([ninf, inf])], &
      [real(dp) :: nan])
   call check('all missing decisive', [all_value([nan, 0.0_dp], value = 1.0_dp)], &
      [real(dp) :: 0.0_dp])
   call check('any missing decisive', [any_value([nan, 1.0_dp], value = 1.0_dp)], &
      [real(dp) :: 1.0_dp])
   call check('count NA', [count_value([nan, 1.0_dp, nan], value = nan)], &
      [real(dp) :: 2.0_dp])
   call check_matrix('colRanks preserveShape', col_ranks(x, ties_method = 'dense', preserve_shape = .true.), &
      [real(dp) :: 2.0_dp, 3.0_dp, 1.0_dp, 4.0_dp, &
      2.0_dp, 3.0_dp, 1.0_dp, 4.0_dp, &
      2.0_dp, 4.0_dp, 1.0_dp, 3.0_dp, &
      2.0_dp, 4.0_dp, 1.0_dp, 3.0_dp, &
      2.0_dp, 3.0_dp, 1.0_dp, 2.0_dp])
   if (failures > 0) error stop 1
   print '(a,i0,a)', 'matrixStats R references: ', checks, ' passed'
contains
   subroutine check(label, actual, expected)
      character(*), intent(in) :: label !! Description of the upstream R reference case.
      real(dp), intent(in) :: actual(:) !! Computed values, flattened in Fortran/R column-major order.
      real(dp), intent(in) :: expected(:) !! Independently generated matrixStats 1.5.0 reference values.
      integer :: i
      logical :: ok
      if (size(actual) /= size(expected)) then
         print *, label, ': shape mismatch ', size(actual), size(expected)
         error stop 1
      end if
      do i = 1, size(actual)
         if (ieee_is_nan(expected(i))) then
            ok = ieee_is_nan(actual(i))
         else if (abs(expected(i)) == inf) then
            ok = actual(i) == expected(i)
         else
            ok = .not. ieee_is_nan(actual(i))
            if (ok) ok = abs(actual(i) - expected(i)) <= 5.0e-11_dp * max(1.0_dp, abs(expected(i)))
         end if
         if (.not. ok) then
            print *, label, ': value ', i, 'actual ', actual(i), ' expected ', expected(i)
            failures = failures + 1
         end if
      end do
      checks = checks + 1
   end subroutine check
   subroutine check_matrix(label, actual, expected)
      character(*), intent(in) :: label !! Description of the matrix-valued reference case.
      real(dp), intent(in) :: actual(:, :) !! Computed matrix in the documented output orientation.
      real(dp), intent(in) :: expected(:) !! Column-major flattened upstream reference matrix.
      call check(label, reshape(actual, [size(actual)]), expected)
   end subroutine check_matrix
end program test_reference
