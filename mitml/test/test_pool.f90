! SPDX-License-Identifier: GPL-2.0-or-later
! Upstream mitml 0.4-5 (2023-03-08), authored by Simon Grund,
! Alexander Robitzsch, and Oliver Luedtke; upstream license GPL (>= 2).
! Modern free-form Fortran translation for Fortran-from-R-packages.
program test_pool
   use mitml, only : MITML_OK, dp, pool_confint, pool_estimates, pooled_estimates
   implicit none
   real(dp) :: lower(2)
   real(dp) :: qhat(2, 4)
   real(dp) :: uhat(2, 2, 4)
   real(dp) :: upper(2)
   type(pooled_estimates) :: result
   integer :: i
   integer :: status

   qhat(:, 1) = [1.0_dp, 2.0_dp]
   qhat(:, 2) = [1.2_dp, 1.8_dp]
   qhat(:, 3) = [0.9_dp, 2.2_dp]
   qhat(:, 4) = [1.1_dp, 2.1_dp]
   do i = 1, 4
      uhat(:, :, i) = reshape([0.038_dp + 0.002_dp * real(i, dp), 0.006_dp, &
         0.006_dp, 0.087_dp + 0.003_dp * real(i, dp)], [2, 2])
   end do

   call pool_estimates(qhat, result, uhat, 120.0_dp)
   call assert_true(result%status == MITML_OK, "pool status")
   call assert_close(result%estimate(1), 1.05_dp, 1.0e-12_dp, "estimate 1")
   call assert_close(result%estimate(2), 2.025_dp, 1.0e-12_dp, "estimate 2")
   call assert_close(result%total(1, 1), 0.0638333333333333_dp, 1.0e-12_dp, "total variance 1")
   call assert_close(result%total(2, 2), 0.130958333333333_dp, 1.0e-12_dp, "total variance 2")
   call assert_close(result%riv(1), 0.484496124031008_dp, 1.0e-12_dp, "RIV 1")
   call assert_close(result%df(1), 20.7981350286222_dp, 1.0e-8_dp, "df 1")
   call assert_close(result%fmi(2), 0.327130246583132_dp, 1.0e-8_dp, "FMI 2")
   call assert_close(result%p_value(1), 4.55344893e-4_dp, 2.0e-8_dp, "p value 1")
   call pool_confint(result, 0.95_dp, lower, upper, status)
   call assert_true(status == MITML_OK, "confidence interval status")
   call assert_true(all(lower < result%estimate) .and. all(upper > result%estimate), "confidence limits bracket estimates")

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

end program test_pool
