! SPDX-License-Identifier: GPL-2.0-or-later
program test_hypothesis
   use moments, only : dp, moments_test_result, agostino_test, anscombe_test, &
      bonett_test, jarque_test, ALTERNATIVE_LESS, ALTERNATIVE_GREATER, &
      MOMENTS_SUCCESS, MOMENTS_INSUFFICIENT_DATA
   implicit none

   real(dp) :: x(20), short_x(4)
   type(moments_test_result) :: result

   x = [-1.8_dp, -1.5_dp, -1.3_dp, -1.1_dp, -0.9_dp, -0.6_dp, -0.2_dp, &
      0.0_dp, 0.1_dp, 0.4_dp, 0.7_dp, 1.0_dp, 1.2_dp, 1.5_dp, 1.9_dp, &
      2.5_dp, 3.2_dp, 4.5_dp, 6.0_dp, 8.0_dp]

   result = agostino_test(x)
   call check(result%status == MOMENTS_SUCCESS, 'agostino status')
   call check_close(result%statistic, 1.1908974463516393_dp, 1.0e-13_dp, 'agostino skew')
   call check_close(result%transformed, 2.3937655712857038_dp, 1.0e-13_dp, 'agostino z')
   call check_close(result%p_value, 0.016676403722485227_dp, 1.0e-14_dp, 'agostino p')

   result = anscombe_test(x)
   call check(result%status == MOMENTS_SUCCESS, 'anscombe status')
   call check_close(result%statistic, 3.7872562478206806_dp, 1.0e-13_dp, 'anscombe kurtosis')
   call check_close(result%transformed, 1.3825715743911153_dp, 1.0e-13_dp, 'anscombe z')
   call check_close(result%p_value, 0.16679627001808148_dp, 1.0e-14_dp, 'anscombe p')

   result = bonett_test(x)
   call check(result%status == MOMENTS_SUCCESS, 'bonett status')
   call check_close(result%statistic, 1.936_dp, 1.0e-13_dp, 'bonett tau')
   call check_close(result%transformed, 0.6981264570314798_dp, 1.0e-13_dp, 'bonett z')
   call check_close(result%p_value, 0.48509811366258676_dp, 1.0e-14_dp, 'bonett p')

   result = jarque_test(x)
   call check(result%status == MOMENTS_SUCCESS, 'jarque status')
   call check_close(result%statistic, 5.243932758866766_dp, 1.0e-13_dp, 'jarque statistic')
   call check_close(result%p_value, 0.07265984543459343_dp, 1.0e-14_dp, 'jarque p')

   result = agostino_test(x, ALTERNATIVE_LESS)
   call check_close(result%p_value, 0.008338201861242614_dp, 1.0e-14_dp, 'agostino less')
   result = agostino_test(x, ALTERNATIVE_GREATER)
   call check_close(result%p_value, 0.9916617981387574_dp, 1.0e-14_dp, 'agostino greater')

   short_x = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp]
   result = agostino_test(short_x)
   call check(result%status == MOMENTS_INSUFFICIENT_DATA, 'agostino short sample')

   print '(a)', 'test_hypothesis: PASS'

contains

   subroutine check(condition, label)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: label
      if (.not. condition) then
         write(*, '(a,1x,a)') 'FAIL:', trim(label)
         error stop 1
      end if
   end subroutine check

   subroutine check_close(actual, expected, tolerance, label)
      real(dp), intent(in) :: actual, expected, tolerance
      character(len=*), intent(in) :: label
      call check(abs(actual - expected) <= tolerance, label)
   end subroutine check_close

end program test_hypothesis
