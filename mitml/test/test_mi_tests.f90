! SPDX-License-Identifier: GPL-2.0-or-later
! Upstream mitml 0.4-5 (2023-03-08), authored by Simon Grund,
! Alexander Robitzsch, and Oliver Luedtke; upstream license GPL (>= 2).
! Modern free-form Fortran translation for Fortran-from-R-packages.
program test_mi_tests
   use mitml, only : MITML_OK, d1_test, d2_test, d3_test, d4_test, dp, mi_test_result
   implicit none
   real(dp) :: d(4)
   real(dp) :: ll_model(4)
   real(dp) :: ll_model_pooled(4)
   real(dp) :: ll_null(4)
   real(dp) :: ll_null_pooled(4)
   real(dp) :: qhat(2, 4)
   real(dp) :: uhat(2, 2, 4)
   type(mi_test_result) :: result
   integer :: i

   qhat(:, 1) = [1.0_dp, 2.0_dp]
   qhat(:, 2) = [1.2_dp, 1.8_dp]
   qhat(:, 3) = [0.9_dp, 2.2_dp]
   qhat(:, 4) = [1.1_dp, 2.1_dp]
   do i = 1, 4
      uhat(:, :, i) = reshape([0.038_dp + 0.002_dp * real(i, dp), 0.006_dp, &
         0.006_dp, 0.087_dp + 0.003_dp * real(i, dp)], [2, 2])
   end do
   d = [63.4118967452301_dp, 63.0697674418605_dp, 63.7440305635148_dp, 65.2788844621514_dp]

   call d1_test(qhat, uhat, result)
   call assert_true(result%status == MITML_OK, "D1 status")
   call assert_close(result%riv, 0.473179184771364_dp, 1.0e-12_dp, "D1 RIV")
   call assert_close(result%f_value, 21.4889833037799_dp, 1.0e-10_dp, "D1 statistic")
   call assert_close(result%df2, 15.6056903257408_dp, 1.0e-9_dp, "D1 denominator df")
   call d1_test(qhat, uhat, result, 120.0_dp)
   call assert_true(result%status == MITML_OK, "finite-df D1 status")
   call assert_close(result%df2, 12.7255136574096_dp, 1.0e-9_dp, "finite-df D1 denominator df")

   call d2_test(d, 2, result)
   call assert_true(result%status == MITML_OK, "D2 status")
   call assert_close(result%riv, 0.00462239140728182_dp, 1.0e-12_dp, "D2 RIV")
   call assert_close(result%f_value, 31.7834528565380_dp, 1.0e-9_dp, "D2 statistic")
   call assert_close(result%df2, 84259.9177632051_dp, 1.0e-4_dp, "D2 denominator df")

   ll_model = [-100.0_dp, -99.5_dp, -100.2_dp, -99.8_dp]
   ll_null = [-104.0_dp, -103.8_dp, -104.3_dp, -103.7_dp]
   ll_model_pooled = [-100.1_dp, -100.0_dp, -100.05_dp, -99.95_dp]
   ll_null_pooled = [-103.6_dp, -103.5_dp, -103.55_dp, -103.45_dp]
   call d3_test(ll_model, ll_null, ll_model_pooled, ll_null_pooled, 5, 3, result)
   call assert_close(result%riv, 0.958333333333332_dp, 1.0e-12_dp, "D3 RIV")
   call assert_close(result%f_value, 1.78723404255319_dp, 1.0e-12_dp, "D3 statistic")
   call assert_close(result%df2, 9.75047258979207_dp, 1.0e-9_dp, "D3 denominator df")

   call d4_test(ll_model, ll_null, -100.05_dp, -103.55_dp, 5, 3, result)
   call assert_close(result%f_value, 1.78723404255319_dp, 1.0e-12_dp, "D4 statistic")
   call assert_close(result%df2, 25.0548204158790_dp, 1.0e-9_dp, "D4 denominator df")
   call d4_test(ll_model, ll_null, -100.05_dp, -103.55_dp, 5, 3, result, robust_riv=.true.)
   call assert_close(result%riv, 0.116666666666665_dp, 1.0e-12_dp, "robust D4 RIV")
   call assert_close(result%f_value, 3.13432835820896_dp, 1.0e-12_dp, "robust D4 statistic")

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

end program test_mi_tests
