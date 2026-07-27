! SPDX-License-Identifier: GPL-2.0-or-later
! Based on BCC1997 0.1.1, Copyright (C) 2017 Haoran Zhang.
program test_strikes
   use bcc1997
   implicit none

   type(bcc_parameters) :: parameters
   type(bcc_result), allocatable :: results(:)
   real(dp), parameter :: strikes(5) = [80.0_dp, 90.0_dp, 100.0_dp, &
      110.0_dp, 120.0_dp]
   real(dp) :: parity
   integer :: i

   parameters = bcc_parameters(kappa_v=1.2_dp, kappa_r=0.35_dp, &
      theta_v=0.035_dp, theta_r=0.025_dp, sigma_v=0.25_dp, &
      sigma_r=0.08_dp, mu_j=-0.04_dp, sigma_j=0.18_dp, rho=-0.5_dp, &
      lambda=0.15_dp, spot=100.0_dp, strike=100.0_dp, &
      variance0=0.04_dp, rate0=0.02_dp, maturity=0.75_dp)

   call bcc_price_strikes(parameters, strikes, results)
   call assert_true(size(results) == size(strikes), 'wrong result count')

   do i = 1, size(results)
      call assert_true(results(i)%converged, 'strike integration failed')
      parity = parameters%spot - strikes(i) * &
         exp(-parameters%rate0 * parameters%maturity)
      call assert_close(results(i)%call - results(i)%put, parity, 3.0e-12_dp)
      call assert_true(results(i)%call >= 0.0_dp, 'negative call value')
      call assert_true(results(i)%put >= 0.0_dp, 'negative put value')
   end do

   do i = 2, size(results)
      call assert_true(results(i)%call < results(i - 1)%call, &
         'call prices are not decreasing in strike')
      call assert_true(results(i)%put > results(i - 1)%put, &
         'put prices are not increasing in strike')
   end do

   print '(a)', 'test_strikes: PASS'

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
end program test_strikes
