! SPDX-License-Identifier: Artistic-2.0
program test_invariants
   use matrixstats
   use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_positive_inf, ieee_is_nan
   implicit none
   real(dp) :: x(3, 4), a(4), nan, inf, tol
   real(dp), allocatable :: z(:, :), v(:), empty(:, :), means(:)
   integer, allocatable :: counts(:), ix(:)
   integer :: n
   nan = ieee_value(0.0_dp, ieee_quiet_nan)
   inf = ieee_value(0.0_dp, ieee_positive_inf)
   tol = 1.0e-11_dp
   n = 0
   x = reshape([1.0_dp, 2.0_dp, 4.0_dp, 2.0_dp, 4.0_dp, 8.0_dp, &
      3.0_dp, 6.0_dp, 12.0_dp, 4.0_dp, 8.0_dp, 16.0_dp], [3, 4])
   a = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp]
   call expect(all(abs(row_sums2(x) - col_sums2(transpose(x))) < tol), 'sum transpose identity')
   call expect(all(abs(row_vars(x) - col_vars(transpose(x))) < tol), 'variance transpose identity')
   call expect(all(abs(row_mads(x) - col_mads(transpose(x))) < tol), 'MAD transpose identity')
   call expect(abs(sum2([1.0e308_dp, 1.0e308_dp, -1.0e308_dp]) / 1.0e308_dp - 1.0_dp) < tol, &
      'sum avoids intermediate overflow')
   call expect(abs(mean2([1.0e308_dp, 1.0e308_dp]) / 1.0e308_dp - 1.0_dp) < tol, 'finite huge mean')
   call expect(abs(sum2([1.0e16_dp, 1.0_dp, -1.0e16_dp]) - 1.0_dp) < tol, 'compensated cancellation')
   z = row_prods_as_matrix()
   call expect(abs(z(1, 1) - 1.0_dp) < tol, 'direct product avoids intermediate overflow')
   call expect(abs(product([1.0e200_dp, 1.0e200_dp, 1.0e-200_dp, 1.0e-200_dp]) - 1.0_dp) < tol, &
      'log product avoids intermediate overflow')
   call expect(log_sum_exp([-inf, -inf]) == -inf, 'log sum zero terms')
   call expect(log_sum_exp([inf, 1.0_dp]) == inf, 'log sum positive infinity')
   call expect(ieee_is_nan(log_sum_exp([inf, nan])), 'missing propagates before infinity')
   call expect(log_sum_exp([inf, nan], na_rm = .true.) == inf, 'missing removed before infinity')
   call expect(abs(log_sum_exp([1000.0_dp, 1000.0_dp]) - 1000.0_dp - log(2.0_dp)) < tol, 'stable large logs')
   v = row_sums2(x, rows = [3, 1, 3], cols = [4, 2])
   call expect(all(v == [24.0_dp, 6.0_dp, 24.0_dp]), 'duplicate ordered selectors')
   v = row_sums2(x, rows = [-2, 0], cols = [-1, -3])
   call expect(all(v == [6.0_dp, 24.0_dp]), 'negative and zero selectors')
   v = row_sums2(x, rows = [1, ms_na_index])
   call expect(v(1) == 10.0_dp .and. ieee_is_nan(v(2)), 'missing row index')
   v = row_sums2(x, rows = [1, ms_na_index], na_rm = .true.)
   call expect(all(v == [10.0_dp, 0.0_dp]), 'remove missing row index')
   v = row_vars(x, center = [1.0_dp])
   call expect(all(ieee_is_nan(v)), 'short center rejected safely')
   allocate(empty(3, 0))
   call expect(all(row_sums2(empty) == 0.0_dp), 'empty row sums')
   call expect(all(row_prods(empty) == 1.0_dp), 'empty row products')
   call expect(all(ieee_is_nan(row_means2(empty))), 'empty row means')
   call expect(all(row_alls(empty) == 1.0_dp), 'vacuously true empty all')
   call expect(all(row_anys(empty) == 0.0_dp), 'empty any false')
   call expect(all(row_mins(empty) == inf), 'empty minimum positive infinity')
   call expect(all(row_maxs(empty) == -inf), 'empty maximum negative infinity')
   call expect(all(row_log_sum_exps(empty) == -inf), 'empty log sum exp')
   z = col_diffs(empty)
   call expect(all(shape(z) == [2, 0]), 'zero columns preserve difference shape')
   z = row_diffs(x, differences = 20)
   call expect(all(shape(z) == [3, 0]), 'oversized difference order empty')
   z = row_diffs(x, differences = 0)
   call expect(all(z == x), 'zeroth difference identity')
   call expect(all(diff2(a, lag = 2) == [2.0_dp, 2.0_dp]), 'lag two difference')
   call expect(abs(var_diff(a, diff = 0) - 5.0_dp / 3.0_dp) < tol, 'variance no differencing')
   call expect(abs(var_diff(a)) < tol, 'linear sequence zero difference variance')
   call expect(ieee_is_nan(mad_diff([1.0_dp, 2.0_dp])), 'singleton differentiated MAD upstream rule')
   call expect(iqr_diff([1.0_dp, 2.0_dp]) == 0.0_dp, 'singleton differentiated IQR')
   call expect(all_value([2.0_dp, -3.0_dp]) == 1.0_dp, 'default all means nonzero')
   call expect(any_value([2.0_dp, 0.0_dp]) == 1.0_dp, 'default any means nonzero')
   call expect(all_value([2.0_dp, -3.0_dp], value = 1.0_dp) == 0.0_dp, 'explicit numeric comparison')
   call expect(ieee_is_nan(all_value([nan, 1.0_dp])), 'all unresolved missing')
   call expect(ieee_is_nan(any_value([nan, 0.0_dp])), 'any unresolved missing')
   call expect(count_value([nan, 1.0_dp], value = nan) == 1.0_dp, 'missing count target')
   call expect(all(sign_tabulate([-inf, -2.0_dp, 0.0_dp, 2.0_dp, inf, nan]) == [2, 1, 2, 1, 1, 1]), &
      'sign bins include infinities in signed totals')
   counts = bin_counts(a, [0.0_dp, 2.0_dp, 4.0_dp, 6.0_dp])
   call expect(sum(counts) == size(a), 'histogram partition identity')
   call bin_means(a * a, a, [0.0_dp, 2.0_dp, 4.0_dp, 6.0_dp], means, counts)
   call expect(abs(sum(means * real(counts, dp)) - sum(a * a)) < tol, 'binned mean reconstructs sum')
   call bin_means([nan, 2.0_dp], [1.0_dp, 1.0_dp], [0.0_dp, 2.0_dp], means, counts, na_rm = .false.)
   call expect(ieee_is_nan(means(1)) .and. counts(1) == 2, 'bin missing propagation and count')
   ix = index_by_row([3, 4])
   call expect(all(ix == [1, 4, 7, 10, 2, 5, 8, 11, 3, 6, 9, 12]), 'row-major to column-major permutation')
   call expect(all(row_collapse(x, [4, 2, 1]) == [4.0_dp, 4.0_dp, 4.0_dp]), 'collapse chooses diagonal-like entries')
   call expect(all(abs(row_weighted_means(x) - row_means2(x)) < tol), 'default equal weights mean')
   call expect(all(ieee_is_nan(row_weighted_means(x, [1.0_dp, -1.0_dp, 1.0_dp, 1.0_dp]))), &
      'negative matrix mean weights rejected')
   call expect(all(abs(row_weighted_vars(x) - row_vars(x)) < tol), 'default equal weight variance')
   z = col_ranks(x, ties_method = 'first', preserve_shape = .true.)
   call expect(all(z(:, 1) == [1.0_dp, 2.0_dp, 3.0_dp]), 'column rank orientation')
   z = t_tx_op_y(x, [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp], '*')
   call expect(all(z(1, :) == a * a), 'row-wise vector recycling')
   z = x_op_y(x, [1.0_dp, 2.0_dp, 3.0_dp], '*')
   call expect(all(z(:, 1) == [1.0_dp, 4.0_dp, 12.0_dp]), 'column-wise vector recycling')
   call expect(abs(mean2([huge(1.0_dp), huge(1.0_dp), -huge(1.0_dp)]) / huge(1.0_dp) - &
      1.0_dp / 3.0_dp) < tol, 'nonfinite residual correction skipped')
   call expect(all(row_weighted_mads(reshape([1.0e308_dp, 1.0e308_dp], [1, 2])) == 0.0_dp), &
      'unweighted matrix weighted MAD uses safe stats fallback')
   call expect(all(ieee_is_nan(row_weighted_medians(x, [1.0_dp, -1.0_dp, 1.0_dp, 1.0_dp]))), &
      'negative matrix median weights rejected')
   call expect(weighted_median([2.0_dp, 0.0_dp, 2.0_dp], [inf, inf, inf]) == 2.0_dp, &
      'deliberate repair of upstream infinite-weight median index bug')
   call expect(abs(weighted_mad([4.0_dp, 1.0_dp, 7.0_dp, 1.0_dp, 9.0_dp, 2.0_dp], &
      idxs = [1, 2, 6]) - 1.4826_dp) < tol, 'unweighted MAD honors vector subset')
   call expect(abs(weighted_var([4.0_dp, 1.0_dp, 7.0_dp, 1.0_dp, 9.0_dp, 2.0_dp], &
      idxs = [1, 2, 6]) - 7.0_dp / 3.0_dp) < tol, 'unweighted variance honors vector subset')
   z = row_avgs_per_col_set(x, reshape([1, 2], [2, 1]), w = reshape([1.0_dp], [1, 1]))
   call expect(all(ieee_is_nan(z)), 'weight matrix shape rejected safely')
   print '(a,i0,a)', 'matrixStats independent invariants: ', n, ' passed'
contains
   subroutine expect(ok, label)
      logical, intent(in) :: ok !! Whether the independent invariant was satisfied.
      character(*), intent(in) :: label !! Human-readable explanation of the tested mathematical property.
      if (.not. ok) then
         print *, 'FAIL: ', label
         error stop 1
      end if
      n = n + 1
   end subroutine expect
   pure function row_prods_as_matrix() result(y)
      real(dp) :: y(1, 1)
      y(:, 1) = row_prods(reshape([1.0e200_dp, 1.0e200_dp, 1.0e-200_dp, 1.0e-200_dp], [1, 4]))
   end function row_prods_as_matrix
end program test_invariants
