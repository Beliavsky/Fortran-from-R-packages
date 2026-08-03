! SPDX-License-Identifier: GPL-2.0-or-later
program test_validation
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
   use jdmbs
   implicit none

   type(jdmbs_control) :: control
   type(jdmbs_result) :: fit
   real(dp) :: spot(2), mu(2), sigma(2), strike(2), matrix(2, 2)

   spot = [100.0_dp, 90.0_dp]
   mu = 0.02_dp
   sigma = 0.2_dp
   strike = 100.0_dp
   matrix = 0.5_dp

   control%day = 0
   call normal_bs(spot, mu, sigma, strike, fit, control)
   call check(fit%status == jdmbs_invalid_argument, 'invalid day')

   control%day = 30
   sigma(1) = -0.1_dp
   call normal_bs(spot, mu, sigma, strike, fit, control)
   call check(fit%status == jdmbs_invalid_argument, 'negative sigma')

   sigma = 0.2_dp
   spot(1) = ieee_value(spot(1), ieee_quiet_nan)
   call normal_bs(spot, mu, sigma, strike, fit, control)
   call check(fit%status == jdmbs_nonfinite_input, 'NaN input')

   spot = [100.0_dp, 90.0_dp]
   matrix(1, 1) = 1.2_dp
   call jdm_new_bs(matrix, spot, mu, sigma, 1.0_dp, strike, fit, control)
   call check(fit%status == jdmbs_invalid_argument, 'invalid transmission')

   call jdm_bs(spot, mu, sigma, -1.0_dp, strike, fit, control)
   call check(fit%status == jdmbs_invalid_argument, 'negative lambda')

   print '(a)', 'test_validation: PASS'

contains

   subroutine check(condition, label)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: label
      if (.not. condition) then
         write (*, '(a,1x,a)') 'FAIL:', trim(label)
         error stop 1
      end if
   end subroutine check

end program test_validation
