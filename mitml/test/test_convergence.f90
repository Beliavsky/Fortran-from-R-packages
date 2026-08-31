! SPDX-License-Identifier: GPL-2.0-or-later
! Upstream mitml 0.4-5 (2023-03-08), authored by Simon Grund,
! Alexander Robitzsch, and Oliver Luedtke; upstream license GPL (>= 2).
! Modern free-form Fortran translation for Fortran-from-R-packages.
program test_convergence
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use mitml, only : MITML_OK, dp, gelman_rubin, moving_average, reduced_acf, sd_proportion
   implicit none
   real(dp) :: ac
   real(dp) :: chain(2, 8)
   real(dp) :: eff(2)
   real(dp) :: ma(8)
   real(dp) :: rhat(2)
   real(dp) :: sdp(2)
   real(dp) :: series(8)
   integer :: status

   chain(1, :) = [1.0_dp, 1.1_dp, 0.9_dp, 1.05_dp, 1.02_dp, 0.98_dp, 1.08_dp, 0.92_dp]
   chain(2, :) = 2.0_dp
   call gelman_rubin(chain, 2, rhat, status)
   call assert_true(status == MITML_OK, "Rhat status")
   call assert_close(rhat(1), 0.901598483214702_dp, 1.0e-12_dp, "Rhat reference")
   call assert_close(rhat(2), 1.0_dp, 0.0_dp, "constant-chain Rhat")

   series = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, 6.0_dp, 7.0_dp, 8.0_dp]
   call reduced_acf(series, 1, ac, status)
   call assert_close(ac, 0.625_dp, 1.0e-13_dp, "lag-one ACF")
   call reduced_acf(series, 2, ac, status, smooth=1, smoothing_sd=0.5_dp)
   call assert_close(ac, 0.278881284710438_dp, 1.0e-12_dp, "smoothed ACF")
   call moving_average(series, 1, ma, status)
   call assert_true(status == MITML_OK, "moving-average status")
   call assert_close(ma(1), 1.0_dp, 1.0e-13_dp, "moving-average first edge")
   call assert_close(ma(2), 2.0_dp, 1.0e-13_dp, "moving-average interior")
   call assert_close(ma(8), 8.0_dp, 1.0e-13_dp, "moving-average last edge")

   call sd_proportion(chain, sdp, eff, status, max_order=2)
   call assert_true(status == MITML_OK, "SDprop status")
   call assert_true(ieee_is_finite(sdp(1)) .and. sdp(1) > 0.0_dp, "SDprop finite positive")
   call assert_close(sdp(2), 0.0_dp, 0.0_dp, "constant SDprop")
   call assert_close(eff(2), 8.0_dp, 0.0_dp, "constant effective sample size")

contains

   subroutine assert_true(condition, message)
      logical, intent(in) :: condition !! Predicate that must be true for the regression test to pass.
      character(len=*), intent(in) :: message !! Failure description printed before terminating the test.
      if (.not. condition) error stop "FAIL: " // message
   end subroutine assert_true

   subroutine assert_close(actual, expected, tolerance, message)
      real(dp), intent(in) :: actual !! Value computed by the translated routine.
      real(dp), intent(in) :: expected !! Independent reference value.
      real(dp), intent(in) :: tolerance !! Maximum allowed absolute difference.
      character(len=*), intent(in) :: message !! Failure description printed before terminating the test.
      if (abs(actual - expected) > tolerance) error stop "FAIL: " // message
   end subroutine assert_close

end program test_convergence
