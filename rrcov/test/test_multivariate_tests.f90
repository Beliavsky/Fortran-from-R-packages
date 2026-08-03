! SPDX-License-Identifier: GPL-3.0-or-later
program test_multivariate_tests
  use rrcov, only : dp, test_result, hotelling_t2_one_sample, &
    hotelling_t2_two_sample, wilks_test, rrcov_success
  implicit none
  real(dp) :: x(50, 2), y(50, 2), combined(100, 2)
  integer :: grouping(100), i
  type(test_result) :: result

  do i = 1, 50
    x(i, 1) = 0.2_dp + 0.4_dp * sin(0.17_dp * real(i, dp))
    x(i, 2) = -0.1_dp + 0.3_dp * cos(0.21_dp * real(i, dp))
    y(i, 1) = 2.0_dp + 0.4_dp * sin(0.19_dp * real(i, dp))
    y(i, 2) = 1.5_dp + 0.3_dp * cos(0.23_dp * real(i, dp))
  end do
  combined(1:50, :) = x
  combined(51:100, :) = y
  grouping(1:50) = 1
  grouping(51:100) = 2

  call hotelling_t2_one_sample(x, [0.0_dp, 0.0_dp], result)
  call assert_true(result%status == rrcov_success, "one-sample status")
  call assert_true(result%p_value < 0.05_dp, "one-sample significance")

  call hotelling_t2_two_sample(x, y, result)
  call assert_true(result%status == rrcov_success, "two-sample status")
  call assert_true(result%p_value < 1.0e-6_dp, "two-sample significance")

  call wilks_test(combined, grouping, result, approximation="Bartlett")
  call assert_true(result%status == rrcov_success, "Wilks status")
  call assert_true(result%lambda < 0.5_dp, "Wilks lambda")
  call assert_true(result%p_value < 1.0e-6_dp, "Wilks significance")

  print '(a)', "test_multivariate_tests: PASS"
contains
  subroutine assert_true(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message
    if (.not. condition) then
      write(*, '(a)') "FAIL: " // trim(message)
      error stop 1
    end if
  end subroutine assert_true
end program test_multivariate_tests
