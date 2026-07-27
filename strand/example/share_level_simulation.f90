! SPDX-License-Identifier: GPL-3.0-only
! Upstream authors: Jeff Enos, David Kane, and strand contributors.
program share_level_simulation
  use strand
  implicit none
  type(optimizer_config) :: config
  type(simulation_config) :: simulation
  type(day_result) :: result
  real(dp) :: price(1), volume(1), zero(1), ratio(1), adv(1), alpha(1, 2)
  logical :: investable(1)
  integer :: internal_shares(1, 2), external_shares(1, 2)

  allocate(config%strategies(2))
  call set_strategy(config%strategies(1), 'long', 100.0_dp, 1.0_dp, 0.0_dp)
  call set_strategy(config%strategies(2), 'short', 100.0_dp, 0.0_dp, 1.0_dp)
  price = 10.0_dp
  volume = 0.0_dp
  zero = 0.0_dp
  ratio = 1.0_dp
  adv = 10000.0_dp
  alpha(1, :) = [1.0_dp, -1.0_dp]
  investable = .true.
  internal_shares = 0
  external_shares = 0
  simulation%fill_rate_pct_volume = 100.0_dp

  result = simulate_day(config, simulation, price, price, volume, zero, zero, ratio, adv, &
    investable, alpha, internal_shares, external_shares)
  if (.not. result%success) error stop trim(result%message)
  print '(a,2(i6,1x))', 'orders: ', result%order_shares(1, :)
  print '(a,2(i6,1x))', 'market fills: ', result%market_fill_shares(1, :)
  print '(a,2(i6,1x))', 'internal transfers: ', result%transfer_fill_shares(1, :)
contains
  subroutine set_strategy(spec, name, capital, long_weight, short_weight)
    type(strategy_spec), intent(out) :: spec
    character(len=*), intent(in) :: name
    real(dp), intent(in) :: capital, long_weight, short_weight
    spec%name = name
    spec%capital = capital
    spec%ideal_long_weight = long_weight
    spec%ideal_short_weight = short_weight
    spec%target_long_weight = long_weight
    spec%target_short_weight = short_weight
    spec%has_target_weights = .true.
    spec%position_limit_pct_adv = 100.0_dp
    spec%position_limit_pct_lmv = 100.0_dp
    spec%position_limit_pct_smv = 100.0_dp
    spec%trading_limit_pct_adv = 100.0_dp
  end subroutine set_strategy
end program share_level_simulation
