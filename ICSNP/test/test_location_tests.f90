! SPDX-License-Identifier: GPL-2.0-or-later
program test_location_tests
   use icsnp, only : dp, icsnp_ok, HotellingsT2, rank_ctest, rank_ctest_groups, &
      rank_ictest, test_result
   use icsnp_special, only : chi_square_survival, f_survival, normal_quantile, &
      chi_square_quantile
   implicit none
   integer, parameter :: n = 60, p = 2
   real(dp) :: x(n, p), y(n, p), mu(p)
   integer :: groups(n), i
   type(test_result) :: result

   do i = 1, n
      x(i, 1) = sin(0.31_dp * real(i, dp))
      x(i, 2) = cos(0.19_dp * real(i, dp)) + 0.2_dp * x(i, 1)
      y(i, 1) = x(i, 1) + 0.35_dp
      y(i, 2) = x(i, 2) - 0.20_dp
      groups(i) = 1 + mod(i - 1, 3)
   end do
   mu = 0.0_dp

   call check(abs(chi_square_survival(3.841458820694124_dp, 1.0_dp) - 0.05_dp) < 2.0e-10_dp, &
      'chi-square survival')
   call check(abs(f_survival(3.0_dp, 2.0_dp, 10.0_dp) - 0.095367431640625_dp) < 2.0e-10_dp, &
      'F survival')
   call check(abs(normal_quantile(0.975_dp) - 1.959963984540054_dp) < 2.0e-7_dp, &
      'normal quantile')
   call check(abs(chi_square_quantile(0.9_dp, 3.0_dp) - 6.251388631170325_dp) < 2.0e-8_dp, &
      'chi-square quantile')

   call HotellingsT2(x, result, mu=mu, distribution='f')
   call check_test(result, 'HotellingsT2 one sample')
   call HotellingsT2(x, result, y=y, mu=mu, distribution='chi')
   call check_test(result, 'HotellingsT2 two sample')
   call rank_ctest(x, result, mu=mu, scores='sign')
   call check_test(result, 'rank_ctest sign')
   call rank_ctest(x, result, y=y, scores='rank')
   call check_test(result, 'rank_ctest two sample')
   call rank_ctest_groups(x, groups, result, scores='normal')
   call check_test(result, 'rank_ctest groups')
   call rank_ictest(x, result, mu=mu, scores='rank')
   call check_test(result, 'rank_ictest approximation')
   call rank_ictest(x, result, mu=mu, scores='sign', method='permutation', n_simu=20, seed=4)
   call check(result%status == icsnp_ok .and. result%replications == 20, 'rank_ictest permutation')
   print '(a)', 'test_location_tests: PASS'
contains
   subroutine check_test(value, label)
      type(test_result), intent(in) :: value
      character(len=*), intent(in) :: label
      call check(value%status == icsnp_ok .and. value%statistic >= 0.0_dp .and. &
         value%p_value >= 0.0_dp .and. value%p_value <= 1.0_dp, label)
   end subroutine check_test

   subroutine check(condition, label)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: label
      if (.not. condition) then
         write(*, '(a)') 'FAIL: ' // trim(label)
         error stop 1
      end if
   end subroutine check
end program test_location_tests
