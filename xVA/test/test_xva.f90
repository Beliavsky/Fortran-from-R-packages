program test_xva
  use xva
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  implicit none

  integer :: failures

  failures = 0
  call test_core(failures)
  call test_supervisory(failures)
  call test_capital(failures)
  call test_simulation(failures)
  call test_calculator(failures)

  if (failures /= 0) then
    print '(a,i0)', "FAILED tests: ", failures
    error stop 1
  end if
  print '(a)', "All xVA tests passed."

contains

  subroutine test_core(failures)
    integer, intent(inout) :: failures
    type(csa_t) :: csa
    type(trade_t) :: trades(2)
    type(regulatory_data_t) :: reg_data
    real(dp), allocatable :: pd(:)
    real(dp), allocatable :: grid(:)
    real(dp) :: mtm(3)
    real(dp) :: spread(3)
    real(dp) :: times(3)
    real(dp) :: exposure(2)
    real(dp) :: discount(3)
    real(dp) :: value

    mtm = [2.0_dp, -1.0_dp, 3.0_dp]
    call assert_close(calc_ngr(mtm), 0.8_dp, 1.0e-13_dp, "calc_ngr", failures)

    spread = [100.0_dp, 110.0_dp, 120.0_dp]
    times = [0.0_dp, 1.0_dp, 2.0_dp]
    call calc_pd(spread, 0.5_dp, times, pd)
    call assert_close(pd(1), 1.0_dp - exp(-0.022_dp), 1.0e-13_dp, &
      "calc_pd first", failures)
    call assert_close(pd(2), exp(-0.022_dp) - exp(-0.048_dp), 1.0e-13_dp, &
      "calc_pd second", failures)

    exposure = [10.0_dp, 20.0_dp]
    discount = [1.0_dp, 0.95_dp, 0.90_dp]
    value = calc_va(exposure, discount, pd, 0.5_dp)
    call assert_close(value, -0.5_dp * (10.0_dp * pd(1) + &
      20.0_dp * 0.95_dp * pd(2)), 1.0e-13_dp, "calc_va", failures)

    csa%remargin_frequency = 90.0_dp
    csa%mpor_days = 10.0_dp
    call generate_time_grid(csa, 1.0_dp, grid)
    call assert_true(size(grid) == 9, "time grid size", failures)
    call assert_close(grid(1), 0.0_dp, 1.0e-13_dp, "time grid first", failures)
    call assert_close(grid(2), 80.0_dp / 360.0_dp, 1.0e-13_dp, &
      "time grid lookback", failures)
    call assert_close(grid(9), 1.0_dp, 1.0e-13_dp, "time grid last", failures)

    call assert_true(is_eligible_currency("EUR"), "eligible currency", failures)
    call assert_true(.not. is_eligible_currency("XYZ"), "ineligible currency", failures)
    call assert_true(is_investment_grade("BBB-"), "investment grade", failures)
    call assert_true(.not. is_investment_grade("BB+"), "non-investment grade", failures)

    call trades(1)%configure_class("IRDSwap")
    call trades(2)%configure_class("IRDSwap")
    trades%notional = [2.0_dp, 1.0_dp]
    trades%ei = [3.0_dp, 6.0_dp]
    value = calc_effective_maturity(trades, times, "SA-CCR")
    call assert_close(value, 4.0_dp, 1.0e-13_dp, "effective maturity", failures)

    reg_data%pd_counterparty = 0.002_dp
    reg_data%lgd = 0.45_dp
    value = calc_default_capital(100.0_dp, reg_data, 2.0_dp)
    call assert_close(value, 0.3361282208199462_dp, 1.0e-12_dp, &
      "default capital", failures)
  end subroutine test_core

  subroutine test_supervisory(failures)
    integer, intent(inout) :: failures
    type(supervisory_cva_data_t) :: defaults
    type(supervisory_cva_data_t) :: loaded
    logical :: found
    real(dp) :: value

    call default_supervisory_cva_data(defaults)
    call load_supervisory_cva_data("data", loaded)
    call assert_close(loaded%ir_eligible_correlation(1, 2), 0.91_dp, &
      1.0e-13_dp, "loaded IR correlation", failures)
    call assert_close(loaded%cs_tenor_correlation(1, 2), 0.90_dp, &
      1.0e-13_dp, "loaded CS correlation", failures)
    value = rating_weight(loaded, " A ", found)
    call assert_true(found, "rating lookup found", failures)
    call assert_close(value, 0.008_dp, 1.0e-13_dp, "rating weight", failures)
    value = sector_risk_weight(loaded, loaded%sectors(1), .true., found)
    call assert_true(found, "sector lookup found", failures)
    call assert_close(value, 0.005_dp, 1.0e-13_dp, "sector weight", failures)
    call assert_close(defaults%ir_risk_weight_eligible(1), &
      loaded%ir_risk_weight_eligible(1), 1.0e-13_dp, &
      "default and loaded supervisory data", failures)
    call assert_true(size(loaded%cs_correlation_by_sector%values, 1) == 10, &
      "loaded CS sector matrix", failures)
    call assert_true(size(loaded%cs_sector_counterparty_correlation%values, 1) == 8, &
      "loaded counterparty sector matrix", failures)
    call assert_true(size(loaded%hedge_counterparty_correlations) == 3, &
      "loaded hedge correlations", failures)
    call assert_true(size(loaded%commodity_risk_weights) == 11, &
      "loaded commodity risk weights", failures)
    call assert_true(size(loaded%equity_risk_weights) == 13, &
      "loaded equity risk weights", failures)
    call assert_close(loaded%equity_risk_weights(12)%risk_weight, 0.15_dp, &
      1.0e-13_dp, "equity risk weight", failures)
  end subroutine test_supervisory

  subroutine test_capital(failures)
    integer, intent(inout) :: failures
    type(trade_t) :: trades(1)
    type(regulatory_data_t) :: reg_data
    type(supervisory_cva_data_t) :: supervisory
    type(cva_capital_result_t) :: result
    type(cva_sensitivity_t) :: sensitivities

    call default_supervisory_cva_data(supervisory)
    call trades(1)%configure_class("IRDSwap")
    trades(1)%currency = "USD"

    reg_data%counterparty_rating = "A"
    reg_data%counterparty_sector = supervisory%sectors(1)
    reg_data%cva_framework = "BA-CVA"
    call calc_cva_capital(trades, 100.0_dp, reg_data, supervisory, 2.0_dp, result)
    call assert_close(result%total_charge, 100.0_dp / 1.4_dp * 0.005_dp * &
      2.0_dp, 1.0e-13_dp, "BA-CVA", failures)

    reg_data%cva_framework = "STD-CVA"
    call calc_cva_capital(trades, 100.0_dp, reg_data, supervisory, 2.0_dp, result)
    call assert_close(result%total_charge, 3.5421003238445254_dp, 2.0e-8_dp, &
      "STD-CVA", failures)

    reg_data%cva_framework = "SA-CVA"
    reg_data%sa_cva_multiplier = 1.2_dp
    allocate(sensitivities%interest_rate_delta(1), &
      sensitivities%interest_rate_tenors(1), &
      sensitivities%credit_spread_delta(1), &
      sensitivities%credit_spread_tenors(1))
    sensitivities%interest_rate_delta = 1000.0_dp
    sensitivities%interest_rate_tenors = 1.0_dp
    sensitivities%credit_spread_delta = 100.0_dp
    sensitivities%credit_spread_tenors = 0.5_dp
    call calc_cva_capital(trades, 100.0_dp, reg_data, supervisory, 2.0_dp, &
      result, sensitivities)
    call assert_close(result%interest_rate_charge, 13.32_dp, 1.0e-12_dp, &
      "SA-CVA IR", failures)
    call assert_close(result%credit_spread_charge, 0.6_dp, 1.0e-12_dp, &
      "SA-CVA CS", failures)
    call assert_close(result%total_charge, 13.92_dp, 1.0e-12_dp, &
      "SA-CVA total", failures)
  end subroutine test_capital

  subroutine test_simulation(failures)
    integer, intent(inout) :: failures
    type(trade_t) :: trades(1)
    type(csa_t) :: csa
    type(simulation_data_t) :: sim_data
    type(exposure_profile_t) :: profile
    real(dp), allocatable :: time_points(:)
    real(dp), allocatable :: spot(:)
    real(dp), allocatable :: discount(:)
    integer :: i

    call trades(1)%configure_class("IRDSwap")
    trades(1)%notional = 1.0_dp
    trades(1)%currency = "USD"
    trades(1)%si = 0.0_dp
    trades(1)%ei = 2.0_dp
    trades(1)%buy_sell = "Buy"
    trades(1)%pay_leg_rate = 0.025_dp

    csa%threshold_counterparty = 0.05_dp
    csa%threshold_processing_organization = 0.05_dp
    csa%minimum_transfer_counterparty = 0.005_dp
    csa%minimum_transfer_processing_organization = 0.005_dp
    csa%remargin_frequency = 90.0_dp
    csa%mpor_days = 10.0_dp
    call generate_time_grid(csa, 2.0_dp, time_points)
    allocate(spot(size(time_points)), discount(size(time_points)))
    spot = 0.02_dp
    discount = exp(-time_points * spot)
    sim_data%num_simulations = 30
    sim_data%mean_reversion_a = 0.03_dp
    sim_data%volatility = 0.01_dp
    sim_data%pfe_percentile = 0.90_dp
    call calc_simulated_exposure(discount, time_points, spot, csa, trades, &
      sim_data, "SA-CCR", profile, 12345)
    call assert_true(size(profile%ee) == size(time_points), &
      "simulation profile size", failures)
    call assert_true(all(ieee_is_finite(profile%ee)), &
      "simulation finite EE", failures)
    call assert_true(all(profile%pfe >= 0.0_dp), "simulation nonnegative PFE", failures)
    do i = 2, size(profile%eee)
      call assert_true(profile%eee(i) >= profile%eee(i - 1), &
        "EEE nondecreasing", failures)
    end do
  end subroutine test_simulation

  subroutine test_calculator(failures)
    integer, intent(inout) :: failures
    type(trade_t) :: trades(2)
    type(csa_t) :: csa
    type(collateral_t) :: collateral
    type(simulation_data_t) :: sim_data
    type(regulatory_data_t) :: reg_data
    type(curve_t) :: credit_cpty
    type(curve_t) :: credit_po
    type(curve_t) :: funding
    type(curve_t) :: spot
    type(xva_result_t) :: result_1
    type(xva_result_t) :: result_2

    call make_example_inputs(trades, csa, collateral, sim_data, reg_data, &
      credit_po, credit_cpty, funding, spot)
    reg_data%cva_framework = "SA-CVA"
    call xva_calculator(trades, csa, collateral, sim_data, reg_data, credit_po, &
      credit_cpty, funding, spot, 0.45_dp, 0.45_dp, .false., result_1, 777)
    call xva_calculator(trades, csa, collateral, sim_data, reg_data, credit_po, &
      credit_cpty, funding, spot, 0.45_dp, 0.45_dp, .false., result_2, 777)
    call assert_true(ieee_is_finite(result_1%cva_simulated), &
      "calculator finite CVA", failures)
    call assert_true(result_1%ead >= 0.0_dp, "calculator nonnegative EAD", failures)
    call assert_close(result_1%cva_simulated, result_2%cva_simulated, &
      1.0e-14_dp, "calculator reproducible CVA", failures)
    call assert_close(result_1%cva_capital%total_charge, &
      result_2%cva_capital%total_charge, 1.0e-12_dp, &
      "calculator reproducible capital", failures)
    call assert_true(result_1%has_saccr_values, "calculator SA-CCR values", failures)
  end subroutine test_calculator

  subroutine make_example_inputs(trades, csa, collateral, sim_data, reg_data, &
      credit_po, credit_cpty, funding, spot)
    type(trade_t), intent(out) :: trades(2)
    type(csa_t), intent(out) :: csa
    type(collateral_t), intent(out) :: collateral
    type(simulation_data_t), intent(out) :: sim_data
    type(regulatory_data_t), intent(out) :: reg_data
    type(curve_t), intent(out) :: credit_po
    type(curve_t), intent(out) :: credit_cpty
    type(curve_t), intent(out) :: funding
    type(curve_t), intent(out) :: spot

    call trades(1)%configure_class("IRDSwap")
    trades(1)%external_id = "ext1"
    trades(1)%notional = 1.0_dp
    trades(1)%mtm = 0.045_dp
    trades(1)%currency = "USD"
    trades(1)%si = 0.0_dp
    trades(1)%ei = 1.5_dp
    trades(1)%buy_sell = "Sell"
    trades(1)%pay_leg_rate = 0.05_dp

    call trades(2)%configure_class("IRDSwap")
    trades(2)%external_id = "ext2"
    trades(2)%notional = 1.0_dp
    trades(2)%mtm = -0.065_dp
    trades(2)%currency = "USD"
    trades(2)%si = 0.0_dp
    trades(2)%ei = 2.0_dp
    trades(2)%buy_sell = "Buy"
    trades(2)%pay_leg_rate = 0.05_dp

    csa%id = "csa_1"
    csa%threshold_counterparty = 0.07_dp
    csa%threshold_processing_organization = 0.10_dp
    csa%initial_margin_counterparty = 0.03_dp
    csa%initial_margin_processing_organization = 0.02_dp
    csa%minimum_transfer_counterparty = 0.007_dp
    csa%minimum_transfer_processing_organization = 0.01_dp
    csa%mpor_days = 10.0_dp
    csa%remargin_frequency = 90.0_dp
    csa%values_type = "Actual"

    collateral%id = "col_1"
    collateral%csa_id = "csa_1"
    collateral%amount = 0.03_dp
    collateral%collateral_type = "VariationMargin"

    sim_data%pfe_percentile = 0.90_dp
    sim_data%num_simulations = 20
    sim_data%mean_reversion_a = 0.01_dp
    sim_data%volatility = 0.01_dp

    reg_data%ccr_framework = "SA-CCR"
    reg_data%sa_ccr_simplified = ""
    reg_data%cva_framework = "SA-CVA"
    reg_data%sa_cva_multiplier = 1.2_dp
    reg_data%counterparty_sector = &
      "Sovereigns including central banks and multilateral development banks"
    reg_data%ignore_default_charge = .true.
    reg_data%pd_counterparty = 0.002_dp
    reg_data%pd_processing_organization = 0.005_dp
    reg_data%pd_funding = 0.001_dp
    reg_data%lgd = 0.45_dp
    reg_data%return_on_capital = 0.03_dp
    reg_data%counterparty_rating = "A"
    reg_data%mva_days = 10.0_dp
    reg_data%mva_percentile = 0.99_dp

    allocate(credit_cpty%tenors(4), credit_cpty%rates(4), &
      credit_po%tenors(4), credit_po%rates(4), funding%tenors(4), &
      funding%rates(4))
    credit_cpty%tenors = [0.5_dp, 1.0_dp, 2.0_dp, 5.0_dp]
    credit_cpty%rates = [5.0_dp, 15.0_dp, 30.0_dp, 80.0_dp]
    credit_po%tenors = credit_cpty%tenors
    credit_po%rates = [6.0_dp, 17.0_dp, 35.0_dp, 90.0_dp]
    funding%tenors = credit_cpty%tenors
    funding%rates = [4.0_dp, 12.0_dp, 25.0_dp, 60.0_dp]

    allocate(spot%tenors(5), spot%rates(5))
    spot%tenors = [0.002777778_dp, 0.5_dp, 1.0_dp, 2.0_dp, 5.0_dp]
    spot%rates = [0.0013_dp, 0.0064_dp, 0.0097_dp, 0.0099_dp, 0.0159_dp]
  end subroutine make_example_inputs

  subroutine assert_close(actual, expected, tolerance, label, failures)
    real(dp), intent(in) :: actual
    real(dp), intent(in) :: expected
    real(dp), intent(in) :: tolerance
    character(len=*), intent(in) :: label
    integer, intent(inout) :: failures

    if (.not. ieee_is_finite(actual) .or. abs(actual - expected) > tolerance) then
      print '(a,2(1x,es24.16))', "FAIL " // trim(label) // ":", actual, expected
      failures = failures + 1
    end if
  end subroutine assert_close

  subroutine assert_true(condition, label, failures)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    integer, intent(inout) :: failures

    if (.not. condition) then
      print '(a)', "FAIL " // trim(label)
      failures = failures + 1
    end if
  end subroutine assert_true

end program test_xva
