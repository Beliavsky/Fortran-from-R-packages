! SPDX-License-Identifier: GPL-2.0-or-later
program test_stats
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
   use jfe
   implicit none

   real(dp), parameter :: r(5) = [0.02_dp, -0.01_dp, 0.03_dp, -0.02_dp, 0.01_dp]
   real(dp) :: value
   integer :: status

   call check_close(return_annualized(r, 12.0_dp), 0.0722303532497861_dp, 1.0e-13_dp, &
      'annualized return')
   call check_close(mean_absolute_deviation(r), 0.0168_dp, 1.0e-14_dp, 'mean absolute deviation')
   call check_close(value_at_risk(r, 0.4_dp), -0.01_dp, 1.0e-14_dp, 'type-1 VaR')
   call check_close(expected_shortfall(r, 0.4_dp), -0.015_dp, 1.0e-14_dp, 'expected shortfall')
   call check_close(downside_deviation(r), 0.01_dp, 1.0e-14_dp, 'full downside deviation')
   call check_close(downside_deviation(r, method=downside_subset), sqrt(0.00025_dp), &
      1.0e-14_dp, 'subset downside deviation')
   call check_close(downside_potential(r), 0.006_dp, 1.0e-14_dp, 'downside potential')
   call check_close(upside_risk(r, statistic=3), 0.012_dp, 1.0e-14_dp, 'upside potential')
   call check_close(omega_sharpe_ratio(r), 1.0_dp, 1.0e-14_dp, 'omega sharpe')
   call check_close(volatility_skewness(r, statistic=volatility_ratio), 2.8_dp, &
      1.0e-13_dp, 'volatility skewness')
   call check_close(sharpe_ratio(r), 0.28934569330224724_dp, 1.0e-13_dp, 'sharpe ratio')

   value = value_at_risk(r, -0.1_dp, status)
   if (.not. ieee_is_nan(value)) error stop 'FAIL: invalid VaR should be NaN'
   if (status /= jfe_invalid_argument) error stop 'FAIL: invalid VaR status'

   print '(a)', 'test_stats: PASS'

contains

   subroutine check_close(actual, expected, tolerance, label)
      real(dp), intent(in) :: actual, expected, tolerance
      character(*), intent(in) :: label
      if (abs(actual - expected) > tolerance) then
         print '(a,2(1x,es24.16))', 'FAIL: '//trim(label), actual, expected
         error stop 1
      end if
   end subroutine check_close


end program test_stats
