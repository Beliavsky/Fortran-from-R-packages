program trading_demo
  use trading, only : dp, angular_distance, carbon_footprint, &
    dynamic_beta, dynamic_beta_result, euro_lottery_combination_count, &
    sample_entropy, trade_t
  implicit none

  type(dynamic_beta_result) :: beta_result
  type(trade_t) :: swap
  real(dp) :: returns(8, 2)
  real(dp) :: distances(2, 2)
  real(dp) :: benchmark(80)
  real(dp) :: fund(80)
  integer :: i

  returns(:, 1) = [0.01_dp, -0.02_dp, 0.015_dp, 0.005_dp, -0.01_dp, &
    0.02_dp, -0.005_dp, 0.012_dp]
  returns(:, 2) = [0.008_dp, -0.018_dp, 0.013_dp, 0.006_dp, -0.009_dp, &
    0.018_dp, -0.004_dp, 0.011_dp]
  call angular_distance(returns, distances)

  do i = 1, size(benchmark)
    benchmark(i) = 0.01_dp * sin(0.31_dp * real(i, dp))
    fund(i) = 0.001_dp + 1.4_dp * benchmark(i)
  end do
  call dynamic_beta(fund, benchmark, beta_result, em_iterations=100)

  call swap%configure_class("IRDSwap")
  swap%notional = 10000.0_dp
  swap%si = 0.0_dp
  swap%ei = 10.0_dp
  swap%buy_sell = "Buy"

  write(*, '(a,f10.6)') "Angular distance: ", distances(1, 2)
  write(*, '(a,f10.6)') "Sample entropy: ", sample_entropy(returns(:, 1), tolerance=1.0_dp)
  write(*, '(a,f10.6)') "Smoothed beta: ", &
    beta_result%smoothed_state(1, size(benchmark))
  write(*, '(a,f12.2)') "Adjusted notional: ", swap%calc_adjusted_notional()
  write(*, '(a,f10.6)') "Carbon footprint: ", carbon_footprint(&
    [100.0_dp, 200.0_dp, 50.0_dp], [1000.0_dp, 5000.0_dp, 6000.0_dp], &
    [20000.0_dp, 10000.0_dp, 30000.0_dp])
  write(*, '(a,i0)') "Euro combinations: ", euro_lottery_combination_count()
end program trading_demo
