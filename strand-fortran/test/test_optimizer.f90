! SPDX-License-Identifier: GPL-3.0-only
! Upstream authors: Jeff Enos, David Kane, and strand contributors.
program test_optimizer
  use strand
  implicit none
  type(optimizer_config) :: cfg
  type(optimization_result) :: result
  type(factor_constraint), allocatable :: fc(:)
  type(category_constraint), allocatable :: cc(:)
  real(dp) :: price(4), volume(4), alpha(4, 1)
  logical :: investable(4)
  integer :: shares(4, 1)

  call configure(cfg)
  price = 10.0_dp
  volume = 100000.0_dp
  investable = .true.
  shares = 0
  alpha(:, 1) = [2.0_dp, 1.0_dp, -1.0_dp, -2.0_dp]

  result = optimize_portfolio(cfg, price, volume, investable, alpha, shares)
  call assert_true(result%success)
  call assert_true(all(result%order_shares(:, 1) == [50, 0, 0, -50]))
  call assert_close(result%objective, 2000.0_dp, 1.0e-9_dp)

  allocate(fc(1))
  allocate(fc(1)%values(4))
  fc(1)%name = 'neutral factor'
  fc(1)%strategy = 1
  fc(1)%values = [1.0_dp, -1.0_dp, 1.0_dp, -1.0_dp]
  fc(1)%lower_weight = 0.0_dp
  fc(1)%upper_weight = 0.0_dp
  result = optimize_portfolio(cfg, price, volume, investable, alpha, shares, factor_constraints=fc)
  call assert_true(result%success)
  call assert_true(all(result%order_shares(:, 1) == [50, 0, -50, 0]))
  call assert_close(dot_product(fc(1)%values, result%trade_nmv(:, 1)), 0.0_dp, 1.0e-9_dp)

  investable(1) = .false.
  result = optimize_portfolio(cfg, price, volume, investable, alpha, shares)
  call assert_true(result%success)
  call assert_true(result%order_shares(1, 1) == 0)
  call assert_true(result%order_shares(2, 1) == 50)

  investable = .true.
  allocate(cc(1))
  allocate(cc(1)%category(4))
  cc(1)%name = 'category neutrality'
  cc(1)%strategy = 1
  cc(1)%category = [1, 2, 1, 2]
  cc(1)%lower_weight = 0.0_dp
  cc(1)%upper_weight = 0.0_dp
  result = optimize_portfolio(cfg, price, volume, investable, alpha, shares, category_constraints=cc)
  call assert_true(result%success)
  call assert_true(all(result%order_shares(:, 1) == [50, 0, -50, 0]))

  shares(:, 1) = [50, 0, 0, -50]
  alpha(:, 1) = [1.0_dp, 2.0_dp, -2.0_dp, -1.0_dp]
  cfg%turnover_limit = 0.0_dp
  result = optimize_portfolio(cfg, price, volume, investable, alpha, shares)
  call assert_true(result%success)
  call assert_true(all(result%order_shares(:, 1) == 0))

  print '(a)', 'test_optimizer: PASS'
contains
  subroutine configure(config)
    type(optimizer_config), intent(out) :: config
    allocate(config%strategies(1))
    config%strategies(1)%name = 'strategy_1'
    config%strategies(1)%capital = 1000.0_dp
    config%strategies(1)%ideal_long_weight = 0.5_dp
    config%strategies(1)%ideal_short_weight = 0.5_dp
    config%strategies(1)%target_long_weight = 0.5_dp
    config%strategies(1)%target_short_weight = 0.5_dp
    config%strategies(1)%has_target_weights = .true.
    config%strategies(1)%position_limit_pct_adv = 100.0_dp
    config%strategies(1)%position_limit_pct_lmv = 100.0_dp
    config%strategies(1)%position_limit_pct_smv = 100.0_dp
    config%strategies(1)%trading_limit_pct_adv = 100.0_dp
  end subroutine configure
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
end program test_optimizer
