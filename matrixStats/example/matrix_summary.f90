! SPDX-License-Identifier: Artistic-2.0
program matrix_summary
   use matrixstats
   use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
   implicit none
   real(dp) :: x(3, 4), missing
   real(dp), allocatable :: q(:, :), ranked(:, :), averages(:, :)
   missing = ieee_value(0.0_dp, ieee_quiet_nan)
   x = reshape([1.0_dp, 2.0_dp, 3.0_dp, 2.0_dp, 4.0_dp, 6.0_dp, &
      3.0_dp, missing, 9.0_dp, 4.0_dp, 8.0_dp, 12.0_dp], [3, 4])
   print '(a,3f10.4)', 'Row means: ', row_means2(x, na_rm = .true.)
   print '(a,4f10.4)', 'Column medians: ', col_medians(x, na_rm = .true.)
   print '(a,3f10.4)', 'Weighted row means: ', &
      row_weighted_means(x, [1.0_dp, 2.0_dp, 1.0_dp, 2.0_dp], na_rm = .true.)
   q = row_quantiles(x, [0.25_dp, 0.5_dp, 0.75_dp], na_rm = .true.)
   print '(a,3f10.4)', 'First-row quartiles: ', q(1, :)
   ranked = col_ranks(x, ties_method = 'average', preserve_shape = .true.)
   print '(a,3f10.4)', 'First-column ranks: ', ranked(:, 1)
   averages = row_avgs_per_col_set(x, reshape([1, 2, 3, 4], [2, 2]), na_rm = .true.)
   print '(a,3f10.4)', 'First-set means: ', averages(:, 1)
   if (abs(sum2(row_sums2(x, na_rm = .true.)) - sum2(col_sums2(x, na_rm = .true.))) > 1.0e-12_dp) error stop 1
end program matrix_summary
