program test_draws
   use posterior_draws
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, ieee_positive_inf, ieee_is_nan, ieee_is_finite
   implicit none
   real(dp) :: x(5, 2, 2), small(4, 1, 2), rows(3, 4), special(3, 3), nan, inf
   real(dp) :: logical_rows(3, 3), left(2, 2, 2), right(1, 2, 2)
   real(dp), allocatable :: a(:, :, :), b(:, :, :), m(:, :), q(:, :), v(:)
   logical, allocatable :: finite(:, :), infinite(:, :), missing(:, :)
   integer :: i, status, checks
   checks = 0
   nan = ieee_value(0.0_dp, ieee_quiet_nan)
   inf = ieee_value(0.0_dp, ieee_positive_inf)
   x = reshape([(real(i, dp), i = 1, 20)], shape(x))
   small = reshape([(real(i, dp), i = 1, 8)], shape(small))

   call draws_to_matrix(x, m)
   call assert(all(shape(m) == [10, 2]), 'draws matrix shape')
   call equal(reshape(m, [20]), [(real(i, dp), i = 1, 20)], 'draw ordering')
   call matrix_to_draws(m, 2, a, status)
   call assert(status == 0 .and. all(a == x), 'matrix array round trip')
   call matrix_to_draws(m, 3, a, status)
   call assert(status /= 0 .and. size(a) == 0, 'incompatible chain count')
   call matrix_to_draws(m, 0, a, status)
   call assert(status /= 0 .and. size(a) == 0, 'zero chain count')
   call merge_draws_chains(x, a)
   call assert(all(shape(a) == [10, 1, 2]), 'merge shape')
   call equal(reshape(a, [20]), [(real(i, dp), i = 1, 20)], 'merge ordering')
   call split_draws_chains(x, a)
   call assert(all(shape(a) == [2, 4, 2]), 'odd public split shape')
   call equal(reshape(a, [16]), [1.0_dp, 2.0_dp, 6.0_dp, 7.0_dp, 3.0_dp, 4.0_dp, 8.0_dp, 9.0_dp, &
      11.0_dp, 12.0_dp, 16.0_dp, 17.0_dp, 13.0_dp, 14.0_dp, 18.0_dp, 19.0_dp], 'odd split drops last')
   call split_draws_chains(x(:1, :, :), a)
   call assert(all(shape(a) == [0, 4, 2]), 'one iteration split is empty')

   call bind_draws_arrays(x, x, 'iteration', a, status)
   call assert(status == 0 .and. all(shape(a) == [10, 2, 2]), 'bind iterations')
   call assert(all(a(:5, :, :) == x) .and. all(a(6:, :, :) == x), 'bind iteration values')
   call bind_draws_arrays(x, x, 'chain', a, status)
   call assert(status == 0 .and. all(shape(a) == [5, 4, 2]), 'bind chains')
   call bind_draws_arrays(x, x, 'variable', a, status)
   call assert(status == 0 .and. all(shape(a) == [5, 2, 4]), 'bind variables')
   call bind_draws_arrays(x, small, 'variable', a, status)
   call assert(status /= 0 .and. size(a) == 0, 'bind mismatched dimensions')
   call bind_draws_arrays(x, x, 'draw', a, status)
   call assert(status /= 0 .and. size(a) == 0, 'draw binding unsupported for arrays as in R')

   call subset_draws_array(x, a, iteration = [4, 2, 2], chain = [2, 1], unique = .false., status = status)
   call assert(status == 0, 'subset status')
   call equal(reshape(a, [12]), [2.0_dp, 2.0_dp, 4.0_dp, 7.0_dp, 7.0_dp, 9.0_dp, &
      12.0_dp, 12.0_dp, 14.0_dp, 17.0_dp, 17.0_dp, 19.0_dp], 'subset sorts with repeats')
   call subset_draws_array(x, a, iteration = [4, 2, 2], variable = [2, 1], status = status)
   call assert(status == 0 .and. size(a, 1) == 2, 'subset unique default')
   call equal(a(:, 1, 1), [12.0_dp, 14.0_dp], 'variable requested ordering')
   call subset_draws_array(x, a, iteration = [2, 5], exclude = .true., status = status)
   call equal(a(:, 1, 1), [1.0_dp, 3.0_dp, 4.0_dp], 'subset complement')
   call subset_draws_array(x, a, draw = [9, 2, 2], unique = .false., status = status)
   call assert(status == 0 .and. size(a, 2) == 1, 'draw subset merges chains')
   call equal(a(:, 1, 1), [2.0_dp, 2.0_dp, 9.0_dp], 'draw subset sorted')
   call subset_draws_array(x, a, draw = [1], chain = [1], status = status)
   call assert(status /= 0 .and. size(a) == 0, 'conflicting subset indices')
   call subset_draws_array(x, a, variable = [0], status = status)
   call assert(status /= 0 .and. size(a) == 0, 'invalid variable index')
   call subset_draws_array(x, a, iteration = [integer ::], status = status)
   call assert(status == 0 .and. size(a) == 0, 'empty subset')

   call thin_draws_array(x, 1.5_dp, a, status)
   call equal(a(:, 1, 1), [1.0_dp, 2.0_dp, 4.0_dp], 'fractional thinning ties to even')
   call thin_draws_array(x, 2.5_dp, a, status)
   call equal(a(:, 1, 1), [1.0_dp, 4.0_dp], 'fractional thinning rounds 3.5 upward')
   call thin_draws_array(x, 6.0_dp, a, status)
   call assert(status /= 0 .and. size(a) == 0, 'invalid thinning factor')
   call thin_draws_array(x, nan, a, status)
   call assert(status /= 0 .and. size(a) == 0, 'NaN thinning factor')
   call order_draws_array(x, [5, 4, 3, 2, 1], [2, 1], a, status)
   call equal(a(:, 1, 1), [10.0_dp, 9.0_dp, 8.0_dp, 7.0_dp, 6.0_dp], 'sort numeric dimension labels')
   call order_draws_array(x, [1], [2, 1], a, status)
   call assert(status /= 0 .and. size(a) == 0, 'order rejects mismatched labels')

   call resample_draws_array(small, [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp], 7, a, status = status)
   call assert(status == 0 .and. all(shape(a) == [7, 1, 2]), 'deterministic resample shape')
   call equal(a(:, 1, 1), [1.0_dp, 2.0_dp, 3.0_dp, 3.0_dp, 4.0_dp, 4.0_dp, 4.0_dp], &
      'deterministic resample R reference')
   call resample_draws_array(small, [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp], 7, a, method = 'stratified', &
      uniforms = [0.703724499791861_dp, 0.708184336544946_dp, 0.749903351534158_dp, 0.468353937380016_dp], &
      status = status)
   call assert(status == 0, 'stratified resample status')
   call equal(a(:, 1, 1), [2.0_dp, 2.0_dp, 3.0_dp, 3.0_dp, 4.0_dp, 4.0_dp, 4.0_dp], &
      'stratified resample supplied R uniforms')
   call resample_draws_array(small, [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp], 0, a, status = status)
   call assert(status == 0 .and. all(shape(a) == [0, 1, 2]), 'zero resamples')
   call resample_draws_array(small, [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], 2, a, status = status)
   call assert(status /= 0 .and. size(a) == 0, 'zero weight resampling is invalid')
   call resample_draws_array(small, [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp], 2, a, method = 'stratified', status = status)
   call assert(status /= 0 .and. size(a) == 0, 'stratified requires uniforms')
   call extract_draws_variable(x, 2, v, status)
   call equal(v, [(real(i, dp), i = 11, 20)], 'extract scalar variable')
   call extract_draws_variable_matrix(x, 2, m, status)
   call assert(all(shape(m) == [5, 2]) .and. all(m == x(:, :, 2)), 'extract variable chain matrix')
   call extract_draws_variable(x, 3, v, status)
   call assert(status /= 0 .and. size(v) == 0, 'extract invalid variable')

   call weight_draws_array(small, [0.0_dp, 1.0_dp, 2.0_dp, 3.0_dp], a, status = status)
   call assert(status == 0 .and. size(a, 3) == 3 .and. all(a(:, :, :2) == small), 'append log weights')
   call extract_draws_weights(a, 3, v, status = status)
   call equal(v, [0.0_dp, 1.0_dp, 2.0_dp, 3.0_dp] / 6.0_dp, 'normalized weights')
   call extract_draws_weights(a, 3, v, normalize = .false., status = status)
   call equal(v, [0.0_dp, 1.0_dp, 2.0_dp, 3.0_dp], 'recover unnormalized weights')
   call weight_draws_array(a, [1000.0_dp, 1000.0_dp, -inf, 1000.0_dp], b, log_scale = .true., &
      weight_variable = 3, status = status)
   call assert(status == 0 .and. size(b, 3) == 3, 'replace with log weights')
   call extract_draws_weights(b, 3, v, status = status)
   call equal(v, [1.0_dp, 1.0_dp, 0.0_dp, 1.0_dp] / 3.0_dp, 'stable huge log-weight normalization')
   call weight_draws_array(small, [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], a, status = status)
   call extract_draws_weights(a, 3, v, status = status)
   call assert(status == 0 .and. all(ieee_is_nan(v)), 'all zero weights normalize to NaN as R')
   call weight_draws_array(small, [-1.0_dp, 1.0_dp, 2.0_dp, 3.0_dp], a, status = status)
   call assert(status /= 0 .and. size(a) == 0, 'negative linear weight invalid')

   rows(1, :) = [1.0_dp, 3.0_dp, 5.0_dp, 7.0_dp]
   rows(2, :) = [2.0_dp, 4.0_dp, 6.0_dp, 8.0_dp]
   rows(3, :) = [-1.0_dp, 0.0_dp, 1.0_dp, 2.0_dp]
   call rvar_summary(rows, 'mean', v, status = status)
   call equal(v, [4.0_dp, 5.0_dp, 0.5_dp], 'rvar within mean')
   call rvar_summary(rows, 'median', v)
   call equal(v, [4.0_dp, 5.0_dp, 0.5_dp], 'rvar within median')
   call rvar_summary(rows, 'sum', v)
   call equal(v, [16.0_dp, 20.0_dp, 2.0_dp], 'rvar within sum')
   call rvar_summary(rows, 'prod', v)
   call equal(v, [105.0_dp, 384.0_dp, 0.0_dp], 'rvar within product')
   call rvar_summary(rows, 'min', v)
   call equal(v, [1.0_dp, 2.0_dp, -1.0_dp], 'rvar within minimum')
   call rvar_summary(rows, 'max', v)
   call equal(v, [7.0_dp, 8.0_dp, 2.0_dp], 'rvar within maximum')
   call rvar_summary(rows, 'var', v)
   call equal(v, [20.0_dp, 20.0_dp, 5.0_dp] / 3.0_dp, 'rvar within variance')
   call rvar_summary(rows, 'sd', v)
   call equal(v, sqrt([20.0_dp, 20.0_dp, 5.0_dp] / 3.0_dp), 'rvar within sd')
   call rvar_summary(rows, 'mad', v)
   call equal(v, [2.9652_dp, 2.9652_dp, 1.4826_dp], 'rvar within mad')
   call rvar_summary(rows, 'mad', v, constant = 1.0_dp)
   call equal(v, [2.0_dp, 2.0_dp, 1.0_dp], 'rvar custom MAD constant')
   call rvar_summary(rows, 'unknown', v, status = status)
   call assert(status /= 0 .and. all(ieee_is_nan(v)), 'invalid summary method')
   call rvar_quantiles(rows, [0.0_dp, 0.25_dp, 0.5_dp, 1.0_dp], m, status = status)
   call assert(status == 0 .and. all(shape(m) == [3, 4]), 'rvar quantile shape')
   call equal(reshape(m, [12]), [1.0_dp, 2.0_dp, -1.0_dp, 2.5_dp, 3.5_dp, -0.25_dp, &
      4.0_dp, 5.0_dp, 0.5_dp, 7.0_dp, 8.0_dp, 2.0_dp], 'rvar quantiles type seven')
   call rvar_ranges(rows, m)
   call equal(m(:, 1), [1.0_dp, 2.0_dp, -1.0_dp], 'range minimum')
   call equal(m(:, 2), [7.0_dp, 8.0_dp, 2.0_dp], 'range maximum')
   call rvar_expectation(rows, v)
   call equal(v, [2.0_dp, 7.0_dp, 12.0_dp, 17.0_dp] / 3.0_dp, 'expectation across draws')
   call rvar_over_draws_summary(transpose(rows), 'var', v)
   call equal(v, [20.0_dp, 20.0_dp, 5.0_dp] / 3.0_dp, 'variance across draws')
   call rvar_over_draws_summary(transpose(rows), 'sd', v)
   call equal(v, sqrt([20.0_dp, 20.0_dp, 5.0_dp] / 3.0_dp), 'sd across draws')
   call rvar_over_draws_summary(transpose(rows), 'mad', v)
   call equal(v, [2.9652_dp, 2.9652_dp, 1.4826_dp], 'mad across draws')

   special(1, :) = [1.0_dp, nan, 3.0_dp]
   special(2, :) = [inf, inf, inf]
   special(3, :) = [-inf, inf, 3.0_dp]
   call rvar_summary(special, 'mean', v)
   call equal(v, [nan, inf, nan], 'mean NaN and infinity behavior')
   call rvar_summary(special, 'median', v)
   call equal(v, [nan, inf, 3.0_dp], 'median NaN and infinity behavior')
   call rvar_summary(special, 'mad', v)
   call equal(v, [nan, nan, inf], 'MAD NaN and infinity behavior')
   call rvar_summary(special, 'var', v)
   call assert(all(ieee_is_nan(v)), 'variance nonfinite behavior')
   call rvar_quantiles(special, [0.0_dp, 0.25_dp, 0.5_dp, 0.75_dp, 1.0_dp], m, status = status)
   call assert(status == 0 .and. all(ieee_is_nan(m(1, :))), 'quantiles propagate missing row')
   call equal(m(2, :), [inf, inf, inf, inf, inf], 'quantile repeated infinity endpoints')
   call equal(m(3, :), [-inf, -inf, 3.0_dp, inf, inf], 'quantile mixed infinite endpoints')
   call rvar_quantiles(special, [0.0_dp, 0.5_dp, 1.0_dp], m, na_rm = .true., status = status)
   call equal(m(1, :), [1.0_dp, 2.0_dp, 3.0_dp], 'quantiles remove missing values')
   call rvar_quantiles(rows, [-0.1_dp], m, status = status)
   call assert(status /= 0 .and. all(ieee_is_nan(m)), 'invalid quantile probability')
   call rvar_summary(rows(:, :0), 'sum', v)
   call equal(v, [0.0_dp, 0.0_dp, 0.0_dp], 'empty sums')
   call rvar_summary(rows(:, :0), 'prod', v)
   call equal(v, [1.0_dp, 1.0_dp, 1.0_dp], 'empty products')
   call rvar_summary(rows(:, :0), 'min', v)
   call equal(v, [inf, inf, inf], 'empty minima')
   call rvar_summary(rows(:, :0), 'max', v)
   call equal(v, [-inf, -inf, -inf], 'empty maxima')
   call rvar_summary(rows(:, :0), 'mean', v)
   call assert(all(ieee_is_nan(v)), 'empty mean')
   call rvar_summary(rows(:, :1), 'var', v)
   call assert(all(ieee_is_nan(v)), 'single-element sample variance')

   logical_rows(1, :) = [0.0_dp, nan, 1.0_dp]
   logical_rows(2, :) = [0.0_dp, nan, 0.0_dp]
   logical_rows(3, :) = [1.0_dp, nan, 1.0_dp]
   call rvar_summary(logical_rows, 'all', v)
   call equal(v, [0.0_dp, 0.0_dp, nan], 'three-valued all')
   call rvar_summary(logical_rows, 'any', v)
   call equal(v, [1.0_dp, nan, 1.0_dp], 'three-valued any')
   call rvar_summary(logical_rows, 'all', v, na_rm = .true.)
   call equal(v, [0.0_dp, 0.0_dp, 1.0_dp], 'all missing removal')
   call rvar_summary(logical_rows, 'any', v, na_rm = .true.)
   call equal(v, [1.0_dp, 0.0_dp, 1.0_dp], 'any missing removal')
   call rvar_summary(rows, 'any', v)
   call equal(v, [1.0_dp, 1.0_dp, 1.0_dp], 'numeric coercion for any')
   call rvar_probability(logical_rows, v, na_rm = .true., status = status)
   call equal(v, [1.0_dp / 3.0_dp, nan, 2.0_dp / 3.0_dp], 'probability across logical draws')
   call rvar_probability(rows, v, status = status)
   call assert(status /= 0 .and. all(ieee_is_nan(v)), 'probability rejects nonlogical numbers')
   call rvar_cdf(reshape([1.0_dp, 2.0_dp, nan, inf], [4, 1]), [-inf, 1.0_dp, 1.5_dp, 2.0_dp, inf, nan], &
      m, status)
   call assert(status == 0, 'empirical CDF status')
   call equal(m(:, 1), [0.0_dp, 1.0_dp / 3.0_dp, 1.0_dp / 3.0_dp, 2.0_dp / 3.0_dp, 1.0_dp, nan], &
      'empirical CDF removes missing draws')
   call rvar_cdf(reshape([nan], [1, 1]), [0.0_dp], m, status)
   call assert(status /= 0 .and. all(ieee_is_nan(m)), 'CDF rejects all missing sample')
   call rvar_special_masks(special, finite, infinite, missing)
   call assert(count(finite) == 3 .and. count(infinite) == 5 .and. count(missing) == 1, 'IEEE masks')

   left(1, :, :) = reshape([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp], [2, 2])
   left(2, :, :) = 2.0_dp * left(1, :, :)
   right(1, :, :) = reshape([1.0_dp, 0.0_dp, 0.0_dp, 2.0_dp], [2, 2])
   call rvar_matmul(left, right, a, status)
   call assert(status == 0 .and. all(shape(a) == [2, 2, 2]), 'rvar matrix multiply broadcast')
   call equal(reshape(a(1, :, :), [4]), [1.0_dp, 2.0_dp, 6.0_dp, 8.0_dp], 'rvar first matrix product')
   call equal(reshape(a(2, :, :), [4]), [2.0_dp, 4.0_dp, 12.0_dp, 16.0_dp], 'rvar second matrix product')
   call rvar_matmul(left, right(:, :1, :), a, status)
   call assert(status /= 0 .and. size(a) == 0, 'incompatible matrix multiply')
   call rvar_constant([10.0_dp, 20.0_dp, 30.0_dp], q)
   call assert(all(shape(q) == [1, 3]), 'numeric constant one draw')
   call rvar_ifelse(logical_rows, q, reshape([-1.0_dp], [1, 1]), m, status)
   call assert(status == 0, 'ifelse broadcast status')
   call equal(m(1, :), [-1.0_dp, nan, 30.0_dp], 'ifelse first draw')
   call equal(m(2, :), [-1.0_dp, nan, -1.0_dp], 'ifelse second draw')
   call equal(m(3, :), [10.0_dp, nan, 30.0_dp], 'ifelse third draw')
   call rvar_ifelse(rows, q, q, m, status)
   call assert(status /= 0 .and. size(m) == 0, 'ifelse rejects nonlogical values')
   print '(a,i0,a)', 'posterior draw/rvar tests: ', checks, ' checks passed'
contains
   subroutine assert(condition, label)
      logical, intent(in) :: condition !! Assertion that must hold for this deterministic test.
      character(len=*), intent(in) :: label !! Name of the tested behavior, printed on failure.
      checks = checks + 1
      if (.not. condition) then
         print '(a)', 'FAILED: ' // label
         error stop 1
      end if
   end subroutine assert

   subroutine equal(actual, expected, label)
      real(dp), intent(in) :: actual(:) !! Computed values, including supported IEEE special values.
      real(dp), intent(in) :: expected(:) !! R-reference values in matching order.
      character(len=*), intent(in) :: label !! Name of the numerical assertion.
      integer :: j
      logical :: ok
      call assert(size(actual) == size(expected), label // ' size')
      ok = .true.
      do j = 1, size(actual)
         if (ieee_is_nan(expected(j))) then
            ok = ok .and. ieee_is_nan(actual(j))
         else if (.not. ieee_is_finite(expected(j))) then
            ok = ok .and. actual(j) == expected(j)
         else
            ok = ok .and. ieee_is_finite(actual(j))
            ok = ok .and. abs(actual(j) - expected(j)) <= 2.0e-12_dp * max(1.0_dp, abs(expected(j)))
         end if
      end do
      if (.not. ok) then
         print '(a,*(es23.15,1x))', 'actual: ', actual
         print '(a,*(es23.15,1x))', 'expect: ', expected
      end if
      call assert(ok, label)
   end subroutine equal
end program test_draws
