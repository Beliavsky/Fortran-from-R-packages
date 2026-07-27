! SPDX-License-Identifier: GPL-2.0-or-later
! Based on BCC1997 0.1.1, Copyright (C) 2017 Haoran Zhang.
program test_pricing
   use bcc1997
   implicit none

   type(bcc_parameters) :: parameters, no_jump, changed_jump
   type(bcc_result) :: result, result_no_jump, result_changed_jump
   complex(dp) :: characteristic1, characteristic2
   real(dp) :: parity

   parameters = bcc_parameters(kappa_v=0.5_dp, kappa_r=0.0_dp, &
      theta_v=0.025_dp, theta_r=0.0_dp, sigma_v=0.09_dp, &
      sigma_r=1.0e-7_dp, mu_j=0.0_dp, sigma_j=1.0e-7_dp, rho=0.1_dp, &
      lambda=0.0_dp, spot=100.0_dp, strike=100.0_dp, &
      variance0=0.04_dp, rate0=0.01_dp, maturity=1.0_dp)
   result = bcc_price(parameters)
   call assert_true(result%converged, 'stochastic-volatility case failed')
   call assert_close(result%call, 8.6001374312277_dp, 1.0e-8_dp)
   call assert_close(result%put, 7.6051208061445_dp, 1.0e-8_dp)

   parameters = bcc_parameters(kappa_v=1.5_dp, kappa_r=0.4_dp, &
      theta_v=0.04_dp, theta_r=0.03_dp, sigma_v=0.3_dp, sigma_r=0.1_dp, &
      mu_j=-0.05_dp, sigma_j=0.2_dp, rho=-0.6_dp, lambda=0.2_dp, &
      spot=100.0_dp, strike=105.0_dp, variance0=0.04_dp, &
      rate0=0.025_dp, maturity=1.25_dp)
   result = bcc_price(parameters)
   call assert_true(result%converged, 'full BCC case failed')
   call assert_close(result%call, 7.5828223156054_dp, 2.0e-8_dp)
   call assert_close(result%put, 9.3523119356215_dp, 2.0e-8_dp)
   call assert_close(result%probability1, 0.5956534890820_dp, 2.0e-9_dp)
   call assert_close(result%probability2, 0.5107869439720_dp, 2.0e-9_dp)

   parity = parameters%spot - parameters%strike * &
      exp(-parameters%rate0 * parameters%maturity)
   call assert_close(result%call - result%put, parity, 2.0e-12_dp)

   characteristic1 = bcc_characteristic_1(0.7_dp, parameters)
   characteristic2 = bcc_characteristic_2(0.7_dp, parameters)
   call assert_close(real(characteristic1, dp), -0.98048111242875_dp, &
      2.0e-12_dp)
   call assert_close(aimag(characteristic1), -0.12874376078536_dp, &
      2.0e-12_dp)
   call assert_close(real(characteristic2, dp), -0.96917734187032_dp, &
      2.0e-12_dp)
   call assert_close(aimag(characteristic2), -0.09311634131882_dp, &
      2.0e-12_dp)

   no_jump = parameters
   no_jump%lambda = 0.0_dp
   changed_jump = no_jump
   changed_jump%mu_j = 0.4_dp
   changed_jump%sigma_j = 0.8_dp
   result_no_jump = bcc_price(no_jump)
   result_changed_jump = bcc_price(changed_jump)
   call assert_close(result_no_jump%call, result_changed_jump%call, 1.0e-13_dp)
   call assert_close(result_no_jump%put, result_changed_jump%put, 1.0e-13_dp)

   print '(a)', 'test_pricing: PASS'

contains

   subroutine assert_close(actual, expected, tolerance)
      real(dp), intent(in) :: actual, expected, tolerance

      if (abs(actual - expected) > tolerance * max(1.0_dp, abs(expected))) then
         print '(a,3(1x,es24.16))', 'mismatch:', actual, expected, &
            abs(actual - expected)
         error stop 1
      end if
   end subroutine assert_close

   subroutine assert_true(condition, message)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: message

      if (.not. condition) then
         print '(a)', trim(message)
         error stop 1
      end if
   end subroutine assert_true
end program test_pricing
