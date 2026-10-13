! SPDX-License-Identifier: GPL-2.0-or-later
program normality
   use nortest, only: dp, normality_result, normality_success, ad_test, cvm_test, lillie_test, pearson_test, sf_test
   implicit none
   real(dp), parameter :: x(20) = [-1.6_dp, -1.1_dp, -0.9_dp, -0.8_dp, -0.6_dp, -0.5_dp, -0.4_dp, -0.3_dp, &
                                  -0.2_dp, -0.1_dp, 0.0_dp, 0.1_dp, 0.2_dp, 0.3_dp, 0.4_dp, 0.5_dp, &
                                  0.7_dp, 0.9_dp, 1.2_dp, 1.8_dp]
   type(normality_result) :: tests(5)
   character(len = 16), parameter :: names(5) = [character(len = 16) :: &
                                                'Anderson-Darling', 'Cramer-von Mises', 'Lilliefors', &
                                                'Pearson', 'Shapiro-Francia']
   integer :: i

   tests = [ad_test(x), cvm_test(x), lillie_test(x), pearson_test(x), sf_test(x)]
   print '(a)', 'Normality tests: statistic and approximate p-value'
   do i = 1, size(tests)
      if (tests(i)%status /= normality_success) error stop 'Invalid sample'
      print '(a16, 2f14.8)', names(i), tests(i)%statistic, tests(i)%p_value
   end do
end program normality
