! SPDX-License-Identifier: GPL-2.0-or-later
! Upstream mitml 0.4-5 (2023-03-08), authored by Simon Grund,
! Alexander Robitzsch, and Oliver Luedtke; upstream license GPL (>= 2).
! Modern free-form Fortran translation for Fortran-from-R-packages.
program test_cluster_r2
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_quiet_nan, ieee_value
   use mitml, only : MITML_OK, cluster_means, dp, intraclass_correlation, multilevel_r2, multilevel_r2_result
   implicit none
   integer :: cluster(6)
   integer :: group(6)
   real(dp) :: cov_z(1, 1)
   real(dp) :: means(6)
   real(dp) :: mu_z(1)
   real(dp) :: random_cov(1, 1)
   real(dp) :: random_cross(1)
   real(dp) :: x(4, 1)
   real(dp) :: values(6)
   type(multilevel_r2_result) :: r2
   integer :: status

   cluster = [1, 1, 2, 2, 1, 1]
   group = [1, 1, 1, 1, 2, 2]
   values = [1.0_dp, 3.0_dp, 2.0_dp, 4.0_dp, 10.0_dp, ieee_value(0.0_dp, ieee_quiet_nan)]
   call cluster_means(values, cluster, means, status, group=group)
   call assert_true(status == MITML_OK, "cluster status")
   call assert_close(means(1), 2.0_dp, 1.0e-13_dp, "cluster mean 1")
   call assert_close(means(3), 3.0_dp, 1.0e-13_dp, "cluster mean 2")
   call assert_close(means(5), 10.0_dp, 1.0e-13_dp, "nested group mean")
   call cluster_means(values, cluster, means, status, adjusted=.true., group=group)
   call assert_close(means(1), 3.0_dp, 1.0e-13_dp, "leave-one-out mean")
   call assert_true(ieee_is_nan(means(5)), "single observed member gives missing leave-one-out mean")

   x(:, 1) = [0.0_dp, 1.0_dp, 2.0_dp, 3.0_dp]
   mu_z = [1.5_dp]
   cov_z(1, 1) = 5.0_dp / 3.0_dp
   random_cross = [0.1_dp]
   random_cov(1, 1) = 0.2_dp
   call multilevel_r2(x, [2.0_dp], mu_z, cov_z, 0.8_dp, random_cross, random_cov, 1.2_dp, r2, 1.0_dp, 2.0_dp)
   call assert_true(r2%status == MITML_OK, "R2 status")
   call assert_close(r2%rb1, 0.4_dp, 1.0e-13_dp, "RB1")
   call assert_close(r2%rb2, 0.2_dp, 1.0e-13_dp, "RB2")
   call assert_close(r2%sb, 1.0_dp - 2.0_dp / 3.0_dp, 1.0e-13_dp, "SB")
   call assert_close(r2%mvp, 6.66666666666667_dp / 9.75_dp, 1.0e-12_dp, "MVP")
   call assert_close(intraclass_correlation(1.0_dp, 3.0_dp), 0.25_dp, 1.0e-13_dp, "ICC")

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

end program test_cluster_r2
