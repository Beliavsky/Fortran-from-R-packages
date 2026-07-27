! SPDX-License-Identifier: GPL-3.0-only
! Upstream authors: Jeff Enos, David Kane, and strand contributors.
program strand_demo
  use strand
  implicit none
  type(optimizer_config) :: config
  type(optimization_result) :: result
  real(dp) :: price(4), dollar_volume(4), alpha(4, 1)
  logical :: investable(4)
  integer :: shares(4, 1)

  allocate(config%strategies(1))
  config%strategies(1)%name = 'demo'
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

  price = [50.0_dp, 25.0_dp, 20.0_dp, 40.0_dp]
  dollar_volume = 10000000.0_dp
  alpha(:, 1) = [2.0_dp, 1.0_dp, -1.0_dp, -2.0_dp]
  investable = .true.
  shares = 0
  result = optimize_portfolio(config, price, dollar_volume, investable, alpha, shares)
  if (.not. result%success) error stop trim(result%message)

  print '(a)', 'strand-fortran portfolio optimization demo'
  print '(a,f12.2)', 'objective exposure: ', result%objective
  print '(a,4(i8,1x))', 'order shares: ', result%order_shares(:, 1)
  print '(a,f12.2)', 'joint gross order value: ', sum(abs(result%order_nmv_joint))
end program strand_demo
