! SPDX-License-Identifier: GPL-3.0-or-later
program test_discriminant
  use rrcov, only : dp, lda_model, qda_model, lda_classic_fit, lda_cov_fit, &
    lda_predict, qda_classic_fit, qda_cov_fit, qda_predict, confusion_matrix, rrcov_success
  implicit none
  real(dp) :: x(80, 2), error_rate
  integer :: grouping(80), status, i
  integer, allocatable :: predicted(:), labels(:), table(:, :)
  type(lda_model) :: lda
  type(qda_model) :: qda

  do i = 1, 40
    grouping(i) = 1
    x(i, 1) = -2.0_dp + 0.5_dp * sin(0.31_dp * real(i, dp))
    x(i, 2) = -1.0_dp + 0.5_dp * cos(0.23_dp * real(i, dp))
  end do
  do i = 41, 80
    grouping(i) = 2
    x(i, 1) = 2.0_dp + 0.5_dp * sin(0.29_dp * real(i, dp))
    x(i, 2) = 1.0_dp + 0.5_dp * cos(0.19_dp * real(i, dp))
  end do
  x(1, :) = [10.0_dp, -10.0_dp]
  x(80, :) = [-10.0_dp, 10.0_dp]

  call lda_classic_fit(x, grouping, lda)
  call assert_true(lda%status == rrcov_success, "classical LDA fit")
  call lda_predict(lda, x, predicted, status=status)
  call confusion_matrix(grouping, predicted, labels, table, error_rate, status)
  call assert_true(error_rate < 0.10_dp, "classical LDA error")

  call lda_cov_fit(x, grouping, lda, "ogk")
  call lda_predict(lda, x, predicted, status=status)
  call confusion_matrix(grouping, predicted, labels, table, error_rate, status)
  call assert_true(error_rate < 0.10_dp, "robust LDA error")

  call qda_classic_fit(x, grouping, qda)
  call qda_predict(qda, x, predicted, status=status)
  call confusion_matrix(grouping, predicted, labels, table, error_rate, status)
  call assert_true(error_rate < 0.10_dp, "classical QDA error")

  call qda_cov_fit(x, grouping, qda, "ogk")
  call qda_predict(qda, x, predicted, status=status)
  call confusion_matrix(grouping, predicted, labels, table, error_rate, status)
  call assert_true(error_rate < 0.10_dp, "robust QDA error")

  print '(a)', "test_discriminant: PASS"
contains
  subroutine assert_true(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message
    if (.not. condition) then
      write(*, '(a)') "FAIL: " // trim(message)
      error stop 1
    end if
  end subroutine assert_true
end program test_discriminant
