! SPDX-License-Identifier: MIT
! Modern Fortran translation of R4GoodPersonalFinances computational code.
program test_random_returns
  use r4good_personal_finances
  implicit none
  type(portfolio_spec) :: portfolio
  real(dp), allocatable :: returns(:, :)
  real(dp) :: means(2), sds(2)
  integer :: st, j, n

  call create_default_portfolio(portfolio)
  portfolio%expected_return = [0.0449_dp, 0.02_dp]
  portfolio%standard_deviation = [0.15_dp, 0.0_dp]
  n = 50000
  call generate_random_returns(portfolio, n, returns, st, 1234_i8)
  call assert_true(st == r4gpf_success, "random return status")
  do j = 1, 2
    means(j) = sum(returns(:,j)) / real(n,dp)
    sds(j) = sqrt(sum((returns(:,j)-means(j))**2) / real(n-1,dp))
  end do
  call assert_close(means(1), portfolio%expected_return(1), 0.003_dp, "random mean")
  call assert_close(sds(1), portfolio%standard_deviation(1), 0.003_dp, "random sd")
  call assert_close(maxval(abs(returns(:,2)-portfolio%expected_return(2))), 0.0_dp, 1.0e-14_dp, "zero volatility asset")
  print *, "test_random_returns: PASS"
contains
  subroutine assert_true(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message
    if (.not. condition) then
      print *, "FAIL: ", trim(message)
      error stop 1
    end if
  end subroutine assert_true
  subroutine assert_close(actual, expected, tolerance, message)
    real(dp), intent(in) :: actual, expected, tolerance
    character(len=*), intent(in) :: message
    if (abs(actual - expected) > tolerance) then
      print *, "FAIL: ", trim(message), actual, expected
      error stop 1
    end if
  end subroutine assert_close
end program test_random_returns
