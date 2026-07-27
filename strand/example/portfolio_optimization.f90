! SPDX-License-Identifier: GPL-3.0-only
! Upstream authors: Jeff Enos, David Kane, and strand contributors.
program portfolio_optimization
  use strand
  implicit none
  type(optimizer_config) :: config
  type(factor_constraint), allocatable :: constraints(:)
  type(optimization_result) :: result
  real(dp) :: price(4), adv(4), alpha(4, 1)
  logical :: investable(4)
  integer :: shares(4, 1)

  allocate(config%strategies(1), constraints(1))
  allocate(constraints(1)%values(4))
  config%strategies(1)%name = 'market_neutral'
  config%strategies(1)%capital = 1000000.0_dp
  config%strategies(1)%ideal_long_weight = 0.5_dp
  config%strategies(1)%ideal_short_weight = 0.5_dp
  config%strategies(1)%target_long_weight = 0.5_dp
  config%strategies(1)%target_short_weight = 0.5_dp
  config%strategies(1)%has_target_weights = .true.
  config%strategies(1)%position_limit_pct_adv = 100.0_dp
  config%strategies(1)%position_limit_pct_lmv = 100.0_dp
  config%strategies(1)%position_limit_pct_smv = 100.0_dp
  config%strategies(1)%trading_limit_pct_adv = 100.0_dp

  constraints(1)%name = 'factor neutral'
  constraints(1)%strategy = 1
  constraints(1)%values = [1.0_dp, -1.0_dp, 1.0_dp, -1.0_dp]
  constraints(1)%lower_weight = 0.0_dp
  constraints(1)%upper_weight = 0.0_dp

  price = 50.0_dp
  adv = 10000000.0_dp
  alpha(:, 1) = [2.0_dp, 1.0_dp, -1.0_dp, -2.0_dp]
  investable = .true.
  shares = 0
  result = optimize_portfolio(config, price, adv, investable, alpha, shares, constraints)
  if (.not. result%success) error stop trim(result%message)
  print '(a,4(i8,1x))', 'factor-neutral order shares: ', result%order_shares(:, 1)
end program portfolio_optimization
