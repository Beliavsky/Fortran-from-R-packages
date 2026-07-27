! SPDX-License-Identifier: MIT
program test_portfolio
  use rtl, only: dp, lp_result, refinery_result, frontier_result
  use rtl, only: simplex_maximize, refinery_lp, efficient_frontier_statistics
  implicit none

  type(lp_result) :: lp
  type(refinery_result) :: refinery
  type(frontier_result) :: frontier
  real(dp) :: a(3, 2), b(3), c(2), yields(2, 2)
  real(dp) :: covariance(2, 2)

  a = reshape([1.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp], shape(a))
  b = [4.0_dp, 2.0_dp, 3.0_dp]
  c = [3.0_dp, 2.0_dp]
  lp = simplex_maximize(c, a, b)
  call assert_true(lp%status%ok)
  call assert_close(lp%objective, 10.0_dp, 1.0e-12_dp)
  call assert_close(lp%solution(1), 2.0_dp, 1.0e-12_dp)
  call assert_close(lp%solution(2), 2.0_dp, 1.0e-12_dp)

  yields = reshape([0.5_dp, 0.2_dp, 0.3_dp, 0.6_dp], shape(yields))
  refinery = refinery_lp([40.0_dp, 35.0_dp], [2.0_dp, 3.0_dp], &
    [100.0_dp, 80.0_dp], yields, [100.0_dp, 120.0_dp])
  call assert_true(refinery%status%ok)
  call assert_true(all(refinery%slate >= 0.0_dp))
  call assert_close(refinery%profit, dot_product(refinery%margin, refinery%slate), 1.0e-10_dp)

  covariance = reshape([0.04_dp, 0.01_dp, 0.01_dp, 0.09_dp], shape(covariance))
  frontier = efficient_frontier_statistics(10, [0.08_dp, 0.12_dp], covariance, seed=44)
  call assert_true(frontier%status%ok)
  call assert_close(frontier%weights(1, 1), 1.0_dp, 0.0_dp)
  call assert_close(frontier%weights(2, 2), 1.0_dp, 0.0_dp)
  call assert_true(all(abs(sum(frontier%weights, dim=2) - 1.0_dp) < 1.0e-14_dp))
  call assert_true(frontier%minimum_risk_index >= 1)
  call assert_true(frontier%maximum_sharpe_index >= 1)

  print '(a)', 'test_portfolio: PASS'

contains

  subroutine assert_close(actual, expected, tolerance)
    real(dp), intent(in) :: actual, expected, tolerance
    if (abs(actual - expected) > tolerance) then
      print '(a,3es24.15)', 'mismatch: ', actual, expected, abs(actual - expected)
      error stop 1
    end if
  end subroutine assert_close

  subroutine assert_true(condition)
    logical, intent(in) :: condition
    if (.not. condition) error stop 1
  end subroutine assert_true

end program test_portfolio
