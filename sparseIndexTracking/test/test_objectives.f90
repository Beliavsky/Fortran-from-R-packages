program test_objectives
   use sparse_index_tracking, only : dp, tracking_objective, measure_ete, measure_dr, &
                                     measure_hete, measure_hdr
   implicit none

   real(dp) :: x(3, 2), index_returns(3), weights(2)
   real(dp) :: c1, lambda, p, penalty, value

   x = 0.0_dp
   weights = [0.5_dp, 0.5_dp]
   index_returns = [-0.1_dp, 0.05_dp, 0.3_dp]
   lambda = 0.25_dp
   p = 0.1_dp
   c1 = log(1.0_dp + 1.0_dp / p)
   penalty = 3.0_dp / c1 * sum(log(1.0_dp + weights / p))

   value = tracking_objective(x, index_returns, weights, lambda, p, c1, measure_ete, 0.1_dp)
   call assert_close(value, (0.01_dp + 0.0025_dp + 0.09_dp) / lambda + penalty, 'ETE objective')

   value = tracking_objective(x, index_returns, weights, lambda, p, c1, measure_dr, 0.1_dp)
   call assert_close(value, (0.0025_dp + 0.09_dp) / lambda + penalty, 'DR objective')

   value = tracking_objective(x, index_returns, weights, lambda, p, c1, measure_hete, 0.1_dp)
   call assert_close(value, (0.01_dp + 0.0025_dp + 0.05_dp) / lambda + penalty, 'HETE objective')

   value = tracking_objective(x, index_returns, weights, lambda, p, c1, measure_hdr, 0.1_dp)
   call assert_close(value, (0.0025_dp + 0.05_dp) / lambda + penalty, 'corrected HDR objective')

   value = tracking_objective(x, index_returns, weights, lambda, p, c1, measure_hdr, 0.1_dp, .true.)
   call assert_close(value, 0.05_dp / lambda + penalty, 'source-compatible HDR objective')

   print '(a)', 'test_objectives: PASS'

contains

   subroutine assert_close(actual, expected, message)
      real(dp), intent(in) :: actual, expected
      character(len=*), intent(in) :: message

      if (abs(actual - expected) > 5.0e-13_dp * max(1.0_dp, abs(expected))) then
         print '(a,2es24.14)', 'FAIL: ' // trim(message), actual, expected
         error stop 1
      end if
   end subroutine assert_close

end program test_objectives
