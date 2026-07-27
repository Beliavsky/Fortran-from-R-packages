! SPDX-License-Identifier: GPL-3.0-only
program test_projection
  use portfoliooptim, only : dp, portfolio_result, projection_result, &
    portfolio_optim_projection, zi_projection, risk_mad
  implicit none
  real(dp) :: returns(4, 2), probabilities(4), aconstr(2, 2), bconstr(2)
  real(dp) :: lower(2), upper(2), benchmark(2)
  real(dp) :: c(2), a(2, 2), b(2), xhat(2), bproj(2, 2)
  type(portfolio_result) :: result
  type(projection_result) :: zi

  returns(1, :) = [-0.10_dp, 0.02_dp]
  returns(2, :) = [0.04_dp, -0.03_dp]
  returns(3, :) = [0.08_dp, 0.05_dp]
  returns(4, :) = [0.02_dp, 0.01_dp]
  probabilities = 0.25_dp
  aconstr(1, :) = 1.0_dp
  aconstr(2, :) = -1.0_dp
  bconstr = [1.0_dp, -1.0_dp]
  lower = 0.0_dp
  upper = 1.0_dp
  benchmark = [0.70_dp, 0.30_dp]

  result = portfolio_optim_projection(returns, probabilities, 0.0_dp, &
    risk_mad, benchmark, 0.95_dp, aconstr, bconstr, lower, upper, &
    500, 1.0e-7_dp)
  call assert_true(result%converged)
  call assert_close(result%theta(1), 0.0638297872340426_dp, 5.0e-5_dp)
  call assert_close(result%risk, 0.0197872340425532_dp, 5.0e-7_dp)
  call assert_close(sum(result%theta), 1.0_dp, 1.0e-7_dp)

  returns = 0.02_dp
  result = portfolio_optim_projection(returns, probabilities, 0.01_dp, &
    risk_mad, benchmark, 0.95_dp, aconstr, bconstr, lower, upper, &
    500, 1.0e-7_dp)
  call assert_true(result%converged)
  call assert_close(result%theta(1), benchmark(1), 1.0e-6_dp)
  call assert_close(result%risk, 0.0_dp, 1.0e-12_dp)

  c = 0.0_dp
  a(1, :) = 1.0_dp
  a(2, :) = -1.0_dp
  b = [1.0_dp, -1.0_dp]
  xhat = benchmark
  bproj = 0.0_dp
  bproj(1, 1) = 1.0_dp
  bproj(2, 2) = 1.0_dp
  zi = zi_projection(c, a, b, xhat, bproj, 1000, 1.0e-8_dp)
  call assert_true(zi%converged)
  call assert_close(sum(zi%x), 1.0_dp, 1.0e-6_dp)
  call assert_true(all(zi%x >= 0.0_dp))
  print '(a)', 'test_projection: PASS'

contains

  subroutine assert_close(actual, expected, tolerance)
    real(dp), intent(in) :: actual, expected, tolerance
    if (abs(actual - expected) > tolerance) then
      print *, 'mismatch:', actual, expected, abs(actual - expected)
      error stop 1
    end if
  end subroutine assert_close

  subroutine assert_true(condition)
    logical, intent(in) :: condition
    if (.not. condition) error stop 1
  end subroutine assert_true

end program test_projection
