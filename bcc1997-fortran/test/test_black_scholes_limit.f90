! SPDX-License-Identifier: GPL-2.0-or-later
! Based on BCC1997 0.1.1, Copyright (C) 2017 Haoran Zhang.
program test_black_scholes_limit
   use bcc1997
   implicit none

   type(bcc_parameters) :: parameters
   type(bcc_result) :: result, wrapper_result
   real(dp) :: call_bs, put_bs

   parameters = bcc_parameters(kappa_v=0.0_dp, kappa_r=0.0_dp, &
      theta_v=0.0_dp, theta_r=0.0_dp, sigma_v=1.0e-7_dp, &
      sigma_r=1.0e-7_dp, mu_j=0.0_dp, sigma_j=1.0e-7_dp, rho=0.0_dp, &
      lambda=0.0_dp, spot=100.0_dp, strike=100.0_dp, &
      variance0=0.04_dp, rate0=0.01_dp, maturity=1.0_dp)

   result = bcc_price(parameters)
   call black_scholes_price(100.0_dp, 100.0_dp, 0.01_dp, 0.2_dp, &
      1.0_dp, call_bs, put_bs)

   call assert_true(result%converged, 'BCC integration did not converge')
   call assert_close(result%call, 8.4333186901096_dp, 5.0e-9_dp)
   call assert_close(result%put, 7.4383020650264_dp, 5.0e-9_dp)
   call assert_close(result%call, call_bs, 5.0e-8_dp)
   call assert_close(result%put, put_bs, 5.0e-8_dp)

   wrapper_result = bcc(0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0e-7_dp, &
      1.0e-7_dp, 0.0_dp, 1.0e-7_dp, 0.0_dp, 0.0_dp, 100.0_dp, &
      100.0_dp, 0.04_dp, 0.01_dp, 1.0_dp)
   call assert_close(wrapper_result%call, result%call, 1.0e-12_dp)
   call assert_close(wrapper_result%put, result%put, 1.0e-12_dp)

   print '(a)', 'test_black_scholes_limit: PASS'

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
end program test_black_scholes_limit
