program hypothesis_tests
   use matrixtests, only : dp, test_result, row_t_welch, col_oneway_equalvar, row_cosinor
   implicit none
   real(dp) :: x(2, 6), y(2, 6), times(6)
   type(test_result), allocatable :: tests(:)
   integer :: i
   x(1, :) = [1.0_dp, 2.0_dp, 4.0_dp, 5.0_dp, 7.0_dp, 8.0_dp]
   x(2, :) = [3.0_dp, 1.0_dp, 2.0_dp, 6.0_dp, 5.0_dp, 9.0_dp]
   y = x + 1.0_dp
   tests = row_t_welch(x, y)
   do i = 1, size(tests)
      print '(a,i0,a,f9.5,a,f9.5)', 'Row ', i, ': t = ', tests(i)%statistic, ', p = ', tests(i)%pvalue
   end do
   tests = col_oneway_equalvar(transpose(x), [1, 1, 1, 2, 2, 2])
   print '(a,2f9.5)', 'ANOVA p-values: ', tests%pvalue
   times = [0.0_dp, 4.0_dp, 8.0_dp, 12.0_dp, 16.0_dp, 20.0_dp]
   tests = row_cosinor(x, times)
   print '(a,2f9.5)', 'Cosinor amplitudes: ', tests%amplitude
end program hypothesis_tests
