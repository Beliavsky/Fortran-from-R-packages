! SPDX-License-Identifier: GPL-2.0-or-later
! Based on BCC1997 0.1.1, Copyright (C) 2017 Haoran Zhang.
program test_validation
   use bcc1997
   implicit none

   type(bcc_parameters) :: parameters
   type(integration_settings) :: settings
   type(bcc_result) :: result
   real(dp) :: call, put

   parameters = bcc_parameters()
   parameters%sigma_v = 0.0_dp
   result = bcc_price(parameters)
   call assert_true(result%status == 1, 'zero sigma_v was accepted')

   parameters = bcc_parameters()
   parameters%rho = 1.01_dp
   result = bcc_price(parameters)
   call assert_true(result%status == 1, 'invalid rho was accepted')

   parameters = bcc_parameters()
   parameters%mu_j = -1.0_dp
   result = bcc_price(parameters)
   call assert_true(result%status == 1, 'invalid jump mean was accepted')

   parameters = bcc_parameters()
   settings = integration_settings()
   settings%panel_width = -1.0_dp
   result = bcc_price(parameters, settings)
   call assert_true(result%status == 1, 'invalid quadrature setup was accepted')

   parameters = bcc_parameters(spot=90.0_dp, strike=100.0_dp, maturity=0.0_dp)
   result = bcc_price(parameters)
   call assert_close(result%call, 0.0_dp, 0.0_dp)
   call assert_close(result%put, 10.0_dp, 0.0_dp)
   call assert_true(result%converged, 'zero-maturity result not complete')

   call black_scholes_price(100.0_dp, 100.0_dp, 0.05_dp, 0.2_dp, &
      1.0_dp, call, put)
   call assert_close(call, 10.4505835721856_dp, 2.0e-13_dp)
   call assert_close(put, 5.57352602225697_dp, 2.0e-13_dp)

   print '(a)', 'test_validation: PASS'

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
end program test_validation
