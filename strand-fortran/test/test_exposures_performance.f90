! SPDX-License-Identifier: GPL-3.0-only
! Upstream authors: Jeff Enos, David Kane, and strand contributors.
program test_exposures_performance
  use strand
  implicit none
  type(exposure_result) :: exposures
  type(performance_stats) :: stats
  real(dp) :: nmv(4, 2), capital(2), factors(4, 2)
  integer :: categories(4, 1)

  nmv(:, 1) = [100.0_dp, -50.0_dp, 25.0_dp, -75.0_dp]
  nmv(:, 2) = [20.0_dp, 30.0_dp, -10.0_dp, -40.0_dp]
  capital = [1000.0_dp, 500.0_dp]
  factors(:, 1) = [1.0_dp, 2.0_dp, -1.0_dp, 0.5_dp]
  factors(:, 2) = [0.0_dp, 1.0_dp, 1.0_dp, -1.0_dp]
  categories(:, 1) = [10, 10, 20, 20]

  exposures = calculate_exposures(nmv, capital, factors, categories)
  call assert_true(all(shape(exposures%factor) == [2, 2]))
  call assert_close(exposures%factor(1, 1), -0.0625_dp, 1.0e-14_dp)
  call assert_close(exposures%factor(2, 1), 0.0500_dp, 1.0e-14_dp)
  call assert_close(exposures%factor(1, 2), 0.1400_dp, 1.0e-14_dp)
  call assert_close(exposures%factor(2, 2), 0.1200_dp, 1.0e-14_dp)
  call assert_true(exposures%category_count(1) == 2)
  call assert_true(all(exposures%category_level(:, 1) == [10, 20]))
  call assert_close(exposures%category(1, 1, 1), 0.0500_dp, 1.0e-14_dp)
  call assert_close(exposures%category(2, 1, 1), -0.0500_dp, 1.0e-14_dp)
  call assert_close(exposures%category(1, 1, 2), 0.1000_dp, 1.0e-14_dp)
  call assert_close(exposures%category(2, 1, 2), -0.1000_dp, 1.0e-14_dp)

  stats = summarize_performance([10.0_dp, -5.0_dp, 20.0_dp], &
    [100.0_dp, 100.0_dp, 200.0_dp], [0.0_dp, 20.0_dp, -10.0_dp], &
    [20.0_dp, 10.0_dp, 40.0_dp])
  call assert_close(stats%total_pnl, 25.0_dp, 1.0e-14_dp)
  call assert_close(stats%total_return, 0.15_dp, 1.0e-14_dp)
  call assert_close(stats%annualized_return, 12.6_dp, 1.0e-12_dp)
  call assert_close(stats%annualized_volatility, 1.374772708486752_dp, 1.0e-12_dp)
  call assert_close(stats%annualized_sharpe, 9.165151389911681_dp, 1.0e-12_dp)
  call assert_close(stats%max_drawdown, -0.05_dp, 1.0e-14_dp)
  call assert_close(stats%average_gmv, 400.0_dp / 3.0_dp, 1.0e-12_dp)
  call assert_close(stats%average_nmv, 10.0_dp / 3.0_dp, 1.0e-12_dp)
  call assert_close(stats%average_turnover, 70.0_dp / 3.0_dp, 1.0e-12_dp)
  call assert_close(stats%holding_period_months, 0.544217687074830_dp, 1.0e-12_dp)

  print '(a)', 'test_exposures_performance: PASS'
contains
  subroutine assert_close(actual, expected, tolerance)
    real(dp), intent(in) :: actual, expected, tolerance
    if (abs(actual - expected) > tolerance) then
      print *, 'mismatch:', actual, expected, abs(actual - expected)
      error stop 1
    end if
  end subroutine assert_close
  subroutine assert_true(value)
    logical, intent(in) :: value
    if (.not. value) error stop 1
  end subroutine assert_true
end program test_exposures_performance
