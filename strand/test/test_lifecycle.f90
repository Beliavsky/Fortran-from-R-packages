! SPDX-License-Identifier: GPL-3.0-only
! Upstream authors: Jeff Enos, David Kane, and strand contributors.
program test_lifecycle
  use strand
  implicit none
  type(optimizer_config) :: cfg
  type(simulation_config) :: simcfg
  type(simulation_result) :: result
  real(dp) :: start_price(2, 1), end_price(2, 1), volume(2, 1)
  real(dp) :: dividend(2, 1), distribution(2, 1), ratio(2, 1), adv(2, 1)
  real(dp) :: alpha(2, 1, 1), delisting_return(2, 1)
  logical :: investable(2, 1), delisting(2, 1)
  integer :: initial_shares(1, 1), financing_days(2)

  allocate(cfg%strategies(1))
  cfg%strategies(1)%name = 'strategy'
  cfg%strategies(1)%capital = 100.0_dp
  cfg%strategies(1)%ideal_long_weight = 1.0_dp
  cfg%strategies(1)%ideal_short_weight = 0.0_dp
  cfg%strategies(1)%target_long_weight = 1.0_dp
  cfg%strategies(1)%target_short_weight = 0.0_dp
  cfg%strategies(1)%has_target_weights = .true.
  cfg%strategies(1)%position_limit_pct_adv = 100.0_dp
  cfg%strategies(1)%position_limit_pct_lmv = 100.0_dp
  cfg%strategies(1)%position_limit_pct_smv = 100.0_dp
  cfg%strategies(1)%trading_limit_pct_adv = 100.0_dp

  start_price = 10.0_dp
  end_price = 10.0_dp
  volume = 1000.0_dp
  dividend = 0.0_dp
  distribution = 0.0_dp
  ratio = 1.0_dp
  adv = 10000.0_dp
  alpha = 1.0_dp
  investable = .true.
  delisting = .false.
  delisting(2, 1) = .true.
  delisting_return = 0.0_dp
  delisting_return(2, 1) = -0.25_dp
  initial_shares = 0
  financing_days = [1, 1]
  simcfg%fill_rate_pct_volume = 100.0_dp

  result = simulate_portfolio(cfg, simcfg, start_price, end_price, volume, dividend, &
    distribution, ratio, adv, investable, alpha, initial_shares, financing_days, &
    delisting=delisting, delisting_return=delisting_return)
  call assert_true(result%success)
  call assert_true(result%day(1)%end_shares(1, 1) == 10)
  call assert_close(result%day(1)%gross_pnl(1, 1), 0.0_dp, 1.0e-14_dp)
  call assert_true(result%day(2)%start_shares(1, 1) == 10)
  call assert_true(result%day(2)%end_shares(1, 1) == 0)
  call assert_close(result%day(2)%gross_pnl(1, 1), -25.0_dp, 1.0e-14_dp)
  call assert_true(result%final_shares(1, 1) == 0)
  call assert_close(result%gross_pnl(2, 1), -25.0_dp, 1.0e-14_dp)

  print '(a)', 'test_lifecycle: PASS'
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
end program test_lifecycle
