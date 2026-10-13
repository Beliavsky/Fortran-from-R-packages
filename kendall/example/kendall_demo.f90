! SPDX-License-Identifier: GPL-2.0-or-later
program kendall_demo
   use kendall_module, only : dp, kendall, mann_kendall, seasonal_mann_kendall
   use kendall_module, only : kendall_result, kendall_exact, kendall_ok
   implicit none
   real(dp) :: x(6), y(6), seasonal(12)
   type(kendall_result) :: result
   integer :: i

   x = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, 6.0_dp]
   y = [2.0_dp, 1.0_dp, 4.0_dp, 3.0_dp, 6.0_dp, 5.0_dp]
   call kendall(x, y, result)
   call show('Kendall association, continuity-corrected normal test', result)
   call mann_kendall(y, result, kendall_exact)
   call show('Mann-Kendall trend, optional exact no-ties test', result)
   seasonal = [(real(i, dp), i = 1, 12)]
   call seasonal_mann_kendall(seasonal, 3, result)
   call show('Seasonal trend, three seasons and four cycles', result)

contains

   subroutine show(label, result)
      character(len = *), intent(in) :: label !! Description printed above this example's computed statistics.
      type(kendall_result), intent(in) :: result !! Completed test result, required successful for this example.

      if (result%status /= kendall_ok) error stop 'Unexpected failure in the Kendall example'
      print '(a)', label
      print '(a,f12.8,a,es13.5)', '  tau = ', result%tau, ', p = ', result%sl
      print '(a,f10.3,a,f10.3,a,f10.3)', '  S = ', result%s, ', D = ', result%d, ', var(S) = ', result%vars
   end subroutine show

end program kendall_demo
