! SPDX-License-Identifier: GPL-2.0-or-later
! Upstream mitml 0.4-5 (2023-03-08), authored by Simon Grund,
! Alexander Robitzsch, and Oliver Luedtke; upstream license GPL (>= 2).
! Modern free-form Fortran translation for Fortran-from-R-packages.
program test_constraints
   use mitml, only : MITML_OK, dp, mi_test_result, test_linear_constraints, test_transformed_constraints
   implicit none
   real(dp) :: a(1, 2)
   real(dp) :: b(1)
   real(dp) :: qhat(2, 4)
   real(dp) :: tq(1, 4)
   real(dp) :: tu(1, 1, 4)
   real(dp) :: uhat(2, 2, 4)
   type(mi_test_result) :: direct
   type(mi_test_result) :: linear
   integer :: i

   qhat(:, 1) = [1.0_dp, 2.0_dp]
   qhat(:, 2) = [1.2_dp, 1.8_dp]
   qhat(:, 3) = [0.9_dp, 2.2_dp]
   qhat(:, 4) = [1.1_dp, 2.1_dp]
   a(1, :) = [1.0_dp, -1.0_dp]
   b = [-1.0_dp]
   do i = 1, 4
      uhat(:, :, i) = reshape([0.04_dp, 0.005_dp, 0.005_dp, 0.09_dp], [2, 2])
      tq(:, i) = matmul(a, qhat(:, i)) - b
      tu(:, :, i) = matmul(a, matmul(uhat(:, :, i), transpose(a)))
   end do
   call test_linear_constraints(qhat, uhat, a, b, 1, linear)
   call test_transformed_constraints(tq, tu, 1, direct)
   call assert_true(linear%status == MITML_OK, "linear constraint status")
   call assert_close(linear%f_value, direct%f_value, 1.0e-13_dp, "linear D1 matches transformed D1")
   call assert_close(linear%riv, direct%riv, 1.0e-13_dp, "linear D1 RIV")
   call test_linear_constraints(qhat, uhat, a, b, 2, linear)
   call test_transformed_constraints(tq, tu, 2, direct)
   call assert_close(linear%f_value, direct%f_value, 1.0e-13_dp, "linear D2 matches transformed D2")

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

end program test_constraints
