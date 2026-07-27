! SPDX-License-Identifier: GPL-3.0-only
! Upstream authors: Jeff Enos, David Kane, and strand contributors.
program test_simulation
  use strand
  implicit none
  type(optimizer_config) :: cfg
  type(simulation_config) :: simcfg
  type(day_result) :: day
  real(dp) :: price(1), volume(1), dividend(1), distribution(1), ratio(1), adv(1), alpha(1, 2)
  logical :: investable(1)
  integer :: int_shares(1, 2), ext_shares(1, 2)

  allocate(cfg%strategies(2))
  call set_spec(cfg%strategies(1), 'long', 100.0_dp, 1.0_dp, 0.0_dp)
  call set_spec(cfg%strategies(2), 'short', 100.0_dp, 0.0_dp, 1.0_dp)
  price = 10.0_dp
  volume = 0.0_dp
  dividend = 0.0_dp
  distribution = 0.0_dp
  ratio = 1.0_dp
  adv = 10000.0_dp
  investable = .true.
  alpha(1, 1) = 1.0_dp
  alpha(1, 2) = -1.0_dp
  int_shares = 0
  ext_shares = 0
  simcfg%fill_rate_pct_volume = 100.0_dp

  day = simulate_day(cfg, simcfg, price, price, volume, dividend, distribution, ratio, adv, &
    investable, alpha, int_shares, ext_shares)
  call assert_true(day%success)
  call assert_true(all(day%order_shares(1, :) == [10, -10]))
  call assert_true(all(day%market_order_shares(1, :) == [0, 0]))
  call assert_true(all(day%transfer_fill_shares(1, :) == [10, -10]))
  call assert_true(all(day%end_shares(1, :) == [10, -10]))
  call assert_true(all(ext_shares == 0))
  call assert_true(all(int_shares(1, :) == [10, -10]))

  call test_partial_fill()
  print '(a)', 'test_simulation: PASS'
contains
  subroutine set_spec(spec, name, capital, long_weight, short_weight)
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
  end subroutine set_spec

  subroutine test_partial_fill()
    type(optimizer_config) :: local_cfg
    type(simulation_config) :: local_sim
    type(day_result) :: local_day
    real(dp) :: p(2), v(2), zero(2), one(2), dollar_vol(2), a(2, 1)
    logical :: inv(2)
    integer :: ints(2, 1), exts(2, 1)

    allocate(local_cfg%strategies(1))
    call set_spec(local_cfg%strategies(1), 'market', 200.0_dp, 0.5_dp, 0.5_dp)
    p = 10.0_dp
    v = 5.0_dp
    zero = 0.0_dp
    one = 1.0_dp
    dollar_vol = 10000.0_dp
    a(:, 1) = [1.0_dp, -1.0_dp]
    inv = .true.
    ints = 0
    exts = 0
    local_sim%fill_rate_pct_volume = 100.0_dp
    local_sim%transaction_cost_pct = 1.0_dp
    local_day = simulate_day(local_cfg, local_sim, p, p, v, zero, zero, one, dollar_vol, &
      inv, a, ints, exts)
    call assert_true(local_day%success)
    call assert_true(all(local_day%order_shares(:, 1) == [10, -10]))
    call assert_true(all(local_day%market_fill_shares(:, 1) == [5, -5]))
    call assert_close(sum(local_day%trade_costs(:, 1)), 1.0_dp, 1.0e-12_dp)
  end subroutine test_partial_fill

  subroutine assert_close(actual, expected, tolerance)
    real(dp), intent(in) :: actual, expected, tolerance
    if (abs(actual - expected) > tolerance) error stop 1
  end subroutine assert_close
  subroutine assert_true(value)
    logical, intent(in) :: value
    if (.not. value) error stop 1
  end subroutine assert_true
end program test_simulation
