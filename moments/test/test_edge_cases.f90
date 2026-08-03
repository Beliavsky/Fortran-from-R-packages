! SPDX-License-Identifier: GPL-2.0-or-later
program test_edge_cases
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_value, ieee_quiet_nan
   use moments, only : dp, moment, skewness, kurtosis, geary, jarque_test, &
      moments_test_result, MOMENTS_DEGENERATE_DATA, MOMENTS_NONFINITE_DATA
   implicit none

   real(dp) :: constant(10), with_nan(5)
   type(moments_test_result) :: result

   constant = 2.0_dp
   call check(ieee_is_nan(skewness(constant)), 'constant skewness')
   call check(ieee_is_nan(kurtosis(constant)), 'constant kurtosis')
   call check(ieee_is_nan(geary(constant)), 'constant geary')
   result = jarque_test(constant)
   call check(result%status == MOMENTS_DEGENERATE_DATA, 'constant jarque status')

   with_nan = [1.0_dp, 2.0_dp, ieee_value(0.0_dp, ieee_quiet_nan), 4.0_dp, 5.0_dp]
   call check(ieee_is_nan(moment(with_nan)), 'nan retained')
   call check(abs(moment(with_nan, na_rm=.true.) - 3.0_dp) < 1.0e-14_dp, 'nan removed')
   result = jarque_test(with_nan)
   call check(result%status == MOMENTS_NONFINITE_DATA, 'jarque nan status')

   print '(a)', 'test_edge_cases: PASS'

contains

   subroutine check(condition, label)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: label
      if (.not. condition) then
         write(*, '(a,1x,a)') 'FAIL:', trim(label)
         error stop 1
      end if
   end subroutine check

end program test_edge_cases
