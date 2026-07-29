! stockAnalyst-fortran
! Copyright (C) 2022 MaheshP Kumar (original R package)
! Copyright (C) 2026 Fortran port contributors
! SPDX-License-Identifier: GPL-3.0-only
program test_stock_analyst
  use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
  use stock_analyst
  implicit none

  integer :: failures
  real(dp), parameter :: tol = 1.0e-11_dp
  real(dp) :: vector3(3)

  failures = 0

  call check('share_value_using_ddm_1yr', &
    share_value_using_ddm_1yr(0.20_dp, 50.0_dp, 1.0_dp, 0.08_dp), 46.48_dp)
  call check('share_value_using_ddm_n_years', &
    share_value_using_ddm_n_years([3.0_dp, 3.15_dp], 40.0_dp, &
      [1.0_dp, 2.0_dp], 2.0_dp, 0.08_dp), 39.77_dp)
  call check('share_value_ggm_constant_growth d1', &
    share_value_ggm_constant_growth(1.1024_dp, 0.101_dp, 0.06_dp, 1), 26.89_dp)
  call check('share_value_ggm_constant_growth d0', &
    share_value_ggm_constant_growth(1.04_dp, 0.101_dp, 0.06_dp, 0), 26.89_dp)
  call check('share_value_preferred_stock', &
    share_value_preferred_stock(1.0_dp, 0.09_dp), 11.11_dp)
  call check('share_value_ggm_negative_growth negative', &
    share_value_ggm_negative_growth(4.25_dp, 0.12_dp, -0.10_dp), 19.32_dp)
  call check('share_value_ggm_negative_growth positive', &
    share_value_ggm_negative_growth(4.25_dp, 0.12_dp, 0.10_dp), 19.32_dp)
  call check('computing_g_using_ggm', &
    computing_g_using_ggm(2.0_dp, 0.122_dp, 40.0_dp), 0.0686_dp)
  call check('justified_leading_pe', &
    justified_leading_pe(0.09_dp, 0.32_dp, 0.07_dp), 16.0_dp)
  call check('justified_trailing_pe', &
    justified_trailing_pe(0.09_dp, 0.32_dp, 0.07_dp), 17.1_dp)
  call check('computing_r_with_ggm', &
    computing_r_with_ggm(2.363_dp, 0.055_dp, 56.60_dp), 0.0967_dp)
  call check('share_val_using_two_stage_ddm', &
    share_val_using_two_stage_ddm(0.14_dp, 0.097_dp, 10.0_dp, 0.15_dp, 0.08_dp), 35.9817_dp)
  call check('share_val_using_three_stage_ddm', &
    share_val_using_three_stage_ddm(1.60_dp, 0.12_dp, 2.0_dp, 5.0_dp, &
      0.14_dp, 0.12_dp, 0.102_dp), 224.3515_dp)
  call check('share_val_using_two_stage_hmodel', &
    share_val_using_two_stage_hmodel(0.14_dp, 0.097_dp, 10.0_dp, 5.0_dp, &
      0.15_dp, 0.08_dp), 11.78_dp)
  call check('share_value_no_current_dividend', &
    share_value_no_current_dividend(1.0_dp, 5.0_dp, 0.05_dp, 0.11_dp), 10.98_dp)
  call check('annualized_hpr', &
    annualized_hpr(0.30_dp, 12.0_dp, 10.0_dp, 3.0_dp), 0.0714_dp)
  call check('firm_value_using_disc_fcff', &
    firm_value_using_disc_fcff([0.4_dp, 0.4_dp, 0.4_dp, 0.4_dp], &
      [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp], 0.12_dp), 1.214940_dp)
  call check('equity_value_given_debt_mv', &
    equity_value_given_debt_mv([0.4_dp, 0.4_dp, 0.4_dp, 0.4_dp], &
      [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp], 0.12_dp, 0.21_dp), 1.004940_dp)
  call check('share_value_given_debt_mv', &
    share_value_given_debt_mv([0.4_dp, 0.4_dp, 0.4_dp, 0.4_dp], &
      [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp], 0.12_dp, 0.21_dp, 0.5_dp), 2.01_dp)
  call check('share_value_using_disc_fcfe', &
    share_value_using_disc_fcfe([0.32_dp, 0.34_dp, 0.36_dp], &
      [1.0_dp, 2.0_dp, 3.0_dp], 0.10_dp, 0.5_dp), 1.68_dp)
  call check('firm_value_constant_g', &
    firm_value_constant_g(1.8_dp, 0.08_dp, 0.12_dp), 48.60_dp)
  call check('equity_value_constant_g', &
    equity_value_constant_g(1.8_dp, 0.08_dp, 0.12_dp, 18.0_dp), 30.60_dp)
  call check('share_val_constant_g', &
    share_val_constant_g(1.8_dp, 0.08_dp, 0.12_dp, 1.5_dp), 32.40_dp)
  call check('share_val_two_stage', &
    share_val_two_stage([1.8_dp, 1.8_dp, 1.8_dp, 1.8_dp], &
      [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp], &
      [0.20_dp, 0.20_dp, 0.20_dp, 0.06_dp], 0.124_dp, 0.5_dp), 1.85_dp)
  call check('share_val_three_stage', &
    share_val_three_stage([2.8_dp, 2.8_dp, 2.8_dp, 2.8_dp], &
      [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp], &
      [0.088_dp, 0.088_dp, 0.074_dp, 0.066_dp], 0.10_dp, 0.5_dp), 1.42_dp)
  call check('share_value_ri', &
    share_value_ri(6.0_dp, [1.4_dp, 1.8_dp, 3.175_dp], 0.10_dp, &
      [1.0_dp, 2.0_dp, 3.0_dp]), 11.15_dp)
  call check('share_value_computed_ri', &
    share_value_computed_ri([6.0_dp, 7.0_dp, 8.25_dp], &
      [2.0_dp, 2.5_dp, 4.0_dp], 0.10_dp, [1.0_dp, 2.0_dp, 3.0_dp]), 11.15_dp)

  vector3 = computing_abs_ri([0.5_dp, 1.5_dp, 2.25_dp], &
    [1.0_dp, 2.2_dp, 2.5_dp], [1.0_dp, 2.2_dp, 2.5_dp], &
    0.12_dp, 0.07_dp, 0.30_dp)
  call check_vector('computing_abs_ri', vector3, [0.1810_dp, 0.6782_dp, 1.1525_dp])

  vector3 = computing_ri([6.0_dp, 7.0_dp, 8.25_dp], &
    [2.0_dp, 2.5_dp, 4.0_dp], 0.10_dp)
  call check_vector('computing_ri', vector3, [1.4_dp, 1.8_dp, 3.175_dp])

  call check('share_value_roe', &
    share_value_roe([0.3333_dp, 0.3571_dp, 0.4848_dp], &
      [6.0_dp, 7.0_dp, 8.25_dp], 0.10_dp, [1.0_dp, 2.0_dp, 3.0_dp]), 11.15_dp)
  call check('single_stage_r', &
    single_stage_r(0.16_dp, 18.81_dp, 0.11_dp, 0.08_dp), 50.16_dp)
  call check('share_value_ri_multi_stage_eps', &
    share_value_ri_multi_stage_eps([6.0_dp, 7.0_dp, 8.25_dp], &
      [2.0_dp, 2.5_dp, 4.0_dp], 0.10_dp, [1.0_dp, 2.0_dp, 3.0_dp], &
      1.1_dp, 3.0_dp), 11.97_dp)
  call check('share_value_ri_multi_stage_roe', &
    share_value_ri_multi_stage_roe([0.3333_dp, 0.3571_dp, 0.4848_dp], &
      [6.0_dp, 7.0_dp, 8.25_dp], 0.10_dp, [1.0_dp, 2.0_dp, 3.0_dp], &
      1.1_dp, 3.0_dp), 11.97_dp)
  call check('share_value_ri_plus_pvtv', &
    share_value_ri_plus_pvtv([6.0_dp, 7.0_dp, 8.25_dp], &
      [2.0_dp, 2.5_dp, 4.0_dp], 0.10_dp, [1.0_dp, 2.0_dp, 3.0_dp], &
      0.6_dp, 3), 15.917_dp)
  call check('trailing_pe_basic_eps', &
    trailing_pe_basic_eps(596.5_dp, 15.1_dp), 39.5_dp)
  call check('trailing_pe_diluted_eps', &
    trailing_pe_diluted_eps(596.5_dp, 15.7_dp), 38.0_dp)
  call check('earning_yield_ep', &
    earning_yield_ep(49.19_dp, 3.14_dp), 0.0638_dp)
  call check('leading_pe_next_4qs', &
    leading_pe_next_4qs(15.0_dp, 0.15_dp, 0.18_dp, 0.18_dp, 0.24_dp), 20.0_dp)
  call check('leading_fy1_pe', &
    leading_fy1_pe(184.15_dp, 16.19_dp), 11.4_dp)
  call check('leading_fy2_pe', &
    leading_fy2_pe(184.15_dp, 18.35_dp), 10.0_dp)
  call check('predicted_pe_on_csr', &
    predicted_pe_on_csr(12.12_dp, 2.25_dp, -0.20_dp, 14.43_dp, &
      0.45_dp, 0.9_dp, 0.08_dp), 14.1_dp)
  call check('forward_peg', forward_peg(43.97_dp, 25.30_dp), 1.74_dp)
  call check('predicted_pe_by_fed_model', &
    predicted_pe_by_fed_model(0.0293_dp), 34.1_dp)
  call check('implied_pe_by_yardeni_model', &
    implied_pe_by_yardeni_model(0.06_dp, 0.2_dp, 0.025_dp, 0.0_dp), 18.2_dp)
  call check('share_price_using_past_pe mean', &
    share_price_using_past_pe('mean', &
      [15.8_dp, 23.1_dp, 10.0_dp, 19.8_dp, 35.8_dp], 203.71_dp), 4258.0_dp)
  call check('share_price_using_past_pe median', &
    share_price_using_past_pe('median', &
      [15.8_dp, 23.1_dp, 10.0_dp, 19.8_dp, 35.8_dp], 203.71_dp), 4033.0_dp)
  call check('pe_for_pass_through_inflation', &
    pe_for_pass_through_inflation(0.03_dp, 0.06_dp, 0.70_dp), 20.8_dp)
  call check('terminal_value_using_pe comparable', &
    terminal_value_using_pe('comparable', 14.3_dp, 3.0_dp, &
      0.45_dp, 0.0715_dp, 0.10_dp), 42.90_dp)
  call check('terminal_value_using_pe ggm', &
    terminal_value_using_pe('GGM', 14.3_dp, 3.0_dp, &
      0.45_dp, 0.0715_dp, 0.10_dp), 50.76_dp)
  call check('computing_bv_per_share', &
    computing_bv_per_share(49000.0_dp, 3396.0_dp, 918.2_dp), 49.67_dp)
  call check('computing_pb trailing', &
    computing_pb('trailing', 49.67_dp, 81.23_dp, 0.12_dp, 0.07_dp, 0.10_dp), 1.64_dp)
  call check('computing_pb ggm', &
    computing_pb('GGM', 49.67_dp, 81.23_dp, 0.12_dp, 0.07_dp, 0.10_dp), 1.67_dp)
  call check('computing_ps trailing', &
    computing_ps('trailing', 20.0_dp, 0.35_dp, 0.9_dp, 10.0_dp, 0.07_dp, 0.09_dp), 2.0_dp)
  call check('computing_ps ggm', &
    computing_ps('GGM', 20.0_dp, 0.35_dp, 0.9_dp, 10.0_dp, 0.07_dp, 0.09_dp), 1.7_dp)
  call check('computing_ev_dollar_val', &
    computing_ev_dollar_val(15008.0_dp, 0.0_dp, 2013.0_dp, 4060.0_dp), 12961.0_dp)
  call check('computing_ev_multiple sales', &
    computing_ev_multiple('sales', 14411.0_dp, 3320.0_dp, 18962.0_dp), 0.76_dp)
  call check('computing_ev_multiple ebitda', &
    computing_ev_multiple('EBITDA', 14411.0_dp, 3320.0_dp, 18962.0_dp), 4.3_dp)
  call check('computing_sustainable_g', &
    computing_sustainable_g(0.60_dp, 0.25_dp), 0.15_dp)
  call check('computing_r_with_capm', &
    computing_r_with_capm(0.049_dp, 0.74_dp, 0.045_dp), 0.0823_dp)
  call check('computing_wacc', &
    computing_wacc(35.0_dp, 65.0_dp, 0.056_dp, 0.127_dp, 0.29_dp), 0.09647_dp)
  call check('computing_r_with_ffm', &
    computing_r_with_ffm(0.041_dp, 1.2_dp, 0.5_dp, 0.8_dp, &
      0.055_dp, 0.02_dp, 0.043_dp), 0.151_dp)
  call check('computing_r_with_hmodel', &
    computing_r_with_hmodel(1.0_dp, 20.0_dp, 10.0_dp, 5.0_dp, 0.10_dp, 0.06_dp), 0.123_dp)

  ! Compatibility names retain the original R spelling.
  call check('compatibility shareValueUsingDDM1yr', &
    shareValueUsingDDM1yr(0.20_dp, 50.0_dp, 1.0_dp, 0.08_dp), 46.48_dp)
  call check('compatibility annulizedHPR', &
    annulizedHPR(0.30_dp, 12.0_dp, 10.0_dp, 3.0_dp), 0.0714_dp)
  call check('compatibility computingRwithCAPM', &
    computingRwithCAPM(0.049_dp, 0.74_dp, 0.045_dp), 0.0823_dp)

  call check('round half even 2.5', round_to(2.5_dp, 0), 2.0_dp)
  call check('round half even 3.5', round_to(3.5_dp, 0), 4.0_dp)
  call check('round half even -2.5', round_to(-2.5_dp, 0), -2.0_dp)

  call check_nan('mismatched discounted vectors', &
    share_value_using_ddm_n_years([1.0_dp, 2.0_dp], 10.0_dp, &
      [1.0_dp], 1.0_dp, 0.1_dp))
  call check_nan('mismatched fcfe growth vectors', &
    share_val_two_stage([1.0_dp, 2.0_dp], [1.0_dp, 2.0_dp], &
      [1.0_dp], 0.1_dp, 1.0_dp))

  if (failures /= 0) then
    write (*, '(a,i0)') 'FAILED tests: ', failures
    error stop 1
  end if
  print '(a)', 'All stockAnalyst tests passed.'

contains

  subroutine check(name, actual, expected)
    character(len=*), intent(in) :: name
    real(dp), intent(in) :: actual, expected
    if (abs(actual - expected) > tol) then
      failures = failures + 1
      write (*, '(a,2(1x,es24.16))') 'FAIL '//trim(name)//':', actual, expected
    end if
  end subroutine check

  subroutine check_vector(name, actual, expected)
    character(len=*), intent(in) :: name
    real(dp), intent(in) :: actual(:), expected(:)
    if (size(actual) /= size(expected)) then
      failures = failures + 1
      write (*, '(a)') 'FAIL '//trim(name)//': size mismatch'
    else if (any(abs(actual - expected) > tol)) then
      failures = failures + 1
      write (*, '(a)') 'FAIL '//trim(name)//': values differ'
      write (*, '(a,*(1x,es16.8))') ' actual:  ', actual
      write (*, '(a,*(1x,es16.8))') ' expected:', expected
    end if
  end subroutine check_vector

  subroutine check_nan(name, actual)
    character(len=*), intent(in) :: name
    real(dp), intent(in) :: actual
    if (.not. ieee_is_nan(actual)) then
      failures = failures + 1
      write (*, '(a,1x,es24.16)') 'FAIL '//trim(name)//': expected NaN, got', actual
    end if
  end subroutine check_nan

end program test_stock_analyst
