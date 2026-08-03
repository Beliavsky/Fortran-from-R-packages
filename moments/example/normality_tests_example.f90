! SPDX-License-Identifier: GPL-2.0-or-later
program normality_tests_example
   use moments, only : dp, moments_test_result, agostino_test, anscombe_test, &
      bonett_test, jarque_test
   implicit none

   real(dp) :: x(20)
   type(moments_test_result) :: result

   x = [-1.8_dp, -1.5_dp, -1.3_dp, -1.1_dp, -0.9_dp, -0.6_dp, -0.2_dp, &
      0.0_dp, 0.1_dp, 0.4_dp, 0.7_dp, 1.0_dp, 1.2_dp, 1.5_dp, 1.9_dp, &
      2.5_dp, 3.2_dp, 4.5_dp, 6.0_dp, 8.0_dp]

   result = agostino_test(x)
   call print_result(result)
   result = anscombe_test(x)
   call print_result(result)
   result = bonett_test(x)
   call print_result(result)
   result = jarque_test(x)
   call print_result(result)

contains

   subroutine print_result(test_result)
      type(moments_test_result), intent(in) :: test_result
      write(*, '(a)') trim(test_result%method)
      write(*, '(a,f12.6,a,f12.6)') '  statistic = ', test_result%statistic, &
         ', p-value = ', test_result%p_value
   end subroutine print_result

end program normality_tests_example
