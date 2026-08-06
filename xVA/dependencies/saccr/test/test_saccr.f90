program test_saccr
  use saccr
  implicit none

  integer :: failures

  failures = 0
  call test_scalar_calculations(failures)
  call test_supervisory_data(failures)
  call test_documented_examples(failures)
  call test_basis_and_volatility(failures)
  call test_replacement_cost(failures)
  call test_methodology(failures)
  call test_csv_calculator(failures)
  call test_modes_and_hedge(failures)
  call test_sold_options(failures)

  if (failures /= 0) then
    write(*, '(a,i0)') "FAILED tests: ", failures
    error stop 1
  end if
  write(*, '(a)') "All SACCR tests passed."

contains

  subroutine test_scalar_calculations(failures)
    integer, intent(inout) :: failures

    call assert_close("CalcEAD", calc_ead(60.0_dp, 500.0_dp), 784.0_dp, &
      1.0e-12_dp, failures)
    call assert_close("CalcPFE negative V-C", &
      calc_pfe(-100.0_dp, 500.0_dp, v=-100.0_dp), &
      452.54162246982315_dp, 1.0e-10_dp, failures)
    call assert_close("CalcPFE simplified", &
      calc_pfe(-100.0_dp, 500.0_dp, simplified=.true.), &
      500.0_dp, 1.0e-12_dp, failures)
    call assert_close("normal factor", calculate_factor_multiplier("USD"), &
      1.0_dp, 0.0_dp, failures)
    call assert_close("basis factor", &
      calculate_factor_multiplier("Basis_CDOR_CORRA_USD"), &
      0.5_dp, 0.0_dp, failures)
    call assert_close("volatility factor", &
      calculate_factor_multiplier("Vol_CDOR_EUR"), &
      5.0_dp, 0.0_dp, failures)
  end subroutine test_scalar_calculations

  subroutine test_supervisory_data(failures)
    integer, intent(inout) :: failures
    type(supervisory_record_t), allocatable :: records(:)
    logical :: found
    real(dp) :: value

    call load_supervisory_data("data/supervisory_factors.csv", records)
    call assert_integer("supervisory row count", size(records), 16, failures)
    value = supervisory_factor(records, "CreditSingle", "BBB", found)
    call assert_true("BBB supervisory record found", found, failures)
    call assert_close("BBB supervisory factor", value, 0.0054_dp, &
      1.0e-14_dp, failures)
    value = supervisory_option_volatility(records, "Commodity", "Electricity", found)
    call assert_true("electricity volatility found", found, failures)
    call assert_close("electricity option volatility", value, 1.5_dp, &
      1.0e-14_dp, failures)
  end subroutine test_supervisory_data

  subroutine test_documented_examples(failures)
    integer, intent(inout) :: failures
    type(portfolio_result_t) :: result

    call example_ird(result)
    call assert_close("IRD example", result%total_ead, &
      569.4701409373457_dp, 1.0e-8_dp, failures)

    call example_credit(result)
    call assert_close("Credit example", result%total_ead, &
      381.2383187481058_dp, 1.0e-8_dp, failures)

    call example_commodity(result)
    call assert_close("Commodity example", result%total_ead, &
      5405.615982460716_dp, 1.0e-8_dp, failures)

    call example_fx(result)
    call assert_close("FX example", result%total_ead, &
      924.0_dp, 1.0e-10_dp, failures)

    call example_ird_credit(result)
    call assert_close("IRD plus credit example", result%total_ead, &
      936.4505055406475_dp, 1.0e-8_dp, failures)

    call example_ird_commodity_margined(result)
    call assert_close("margined IRD plus commodity example", result%total_ead, &
      1879.212631496718_dp, 1.0e-8_dp, failures)
  end subroutine test_documented_examples

  subroutine test_basis_and_volatility(failures)
    integer, intent(inout) :: failures
    type(portfolio_result_t) :: result
    logical :: found_basis
    logical :: found_vol
    integer :: i

    call example_basis_volatility(result)
    call assert_close("basis and volatility example", result%total_ead, &
      3522.264264019598_dp, 1.0e-8_dp, failures)
    found_basis = .false.
    found_vol = .false.
    do i = 1, size(result%exposures(1)%addon%hedging_sets)
      if (starts_with(result%exposures(1)%addon%hedging_sets(i)%name, "Basis_")) then
        found_basis = .true.
      end if
      if (starts_with(result%exposures(1)%addon%hedging_sets(i)%name, "Vol_")) then
        found_vol = .true.
      end if
    end do
    call assert_true("basis hedging set present", found_basis, failures)
    call assert_true("volatility hedging set present", found_vol, failures)
  end subroutine test_basis_and_volatility

  subroutine test_replacement_cost(failures)
    integer, intent(inout) :: failures
    type(trade_t) :: trades(2)
    type(csa_t) :: agreement
    type(collateral_t) :: collateral(2)
    type(replacement_cost_t) :: result

    call trades(1)%configure_class("IRDSwap")
    trades(1)%mtm = 100.0_dp
    call trades(2)%configure_class("IRDSwap")
    trades(2)%mtm = -40.0_dp

    call calc_rc(trades, result)
    call assert_close("unmargined V", result%v, 60.0_dp, 0.0_dp, failures)
    call assert_close("unmargined RC", result%rc, 60.0_dp, 0.0_dp, failures)

    agreement%id = "csa_1"
    agreement%minimum_transfer_counterparty = 5.0_dp
    agreement%values_type = "Actual"
    collateral(1)%csa_id = "csa_1"
    collateral(1)%amount = 150.0_dp
    collateral(1)%collateral_type = "ICA"
    collateral(2)%csa_id = "csa_1"
    collateral(2)%amount = 50.0_dp
    collateral(2)%collateral_type = "VariationMargin"
    call calc_rc(trades, result, agreement, collateral)
    call assert_close("collateralized V-C", result%v_c, -140.0_dp, &
      0.0_dp, failures)
    call assert_close("collateralized RC", result%rc, 0.0_dp, &
      0.0_dp, failures)
  end subroutine test_replacement_cost

  subroutine test_methodology(failures)
    integer, intent(inout) :: failures
    type(trade_t) :: trades(2)
    character(len=str_len) :: methodology

    call trades(1)%configure_class("IRDSwap")
    trades(1)%mtm = 20.0_dp
    call trades(2)%configure_class("Commodity")
    trades(2)%mtm = -40.0_dp
    methodology = determine_ccr_methodology(trades, 10000.0_dp)
    call assert_character("methodology", trim(methodology), &
      "Original Exposure Method", failures)
  end subroutine test_methodology

  subroutine test_csv_calculator(failures)
    integer, intent(inout) :: failures
    type(portfolio_result_t) :: result
    integer :: unit

    open(newunit=unit, file="build/test_fx.csv", status="replace", action="write")
    write(unit, '(a)') "TradeObj,Notional,MtM,ccyPair,Si,Ei,BuySell,external_id"
    write(unit, '(a)') "FxForward,250,0,EUR/USD,0,10,Buy,fx_1"
    close(unit)

    call saccr_calculator("build/test_fx.csv", result)
    call assert_close("CSV calculator", result%total_ead, 14.0_dp, &
      1.0e-10_dp, failures)
  end subroutine test_csv_calculator


  subroutine test_modes_and_hedge(failures)
    integer, intent(inout) :: failures
    type(trade_t) :: trades(1)
    type(exposure_result_t) :: exposure
    type(replacement_cost_t) :: rc
    type(csa_t) :: agreement
    type(fx_hedge_result_t) :: hedge

    call trades(1)%configure_class("FxForward")
    trades(1)%notional = 250.0_dp
    trades(1)%mtm = 0.0_dp
    trades(1)%ccy_pair = "EUR/USD"
    trades(1)%si = 0.0_dp
    trades(1)%ei = 10.0_dp
    trades(1)%buy_sell = "Buy"
    call calculate_exposure(trades, exposure, oem=.true.)
    call assert_close("OEM FX EAD", exposure%ead, 140.0_dp, &
      1.0e-10_dp, failures)

    call apply_fx_hedge(14.0_dp, 10.0_dp, 0.2_dp, 0.3_dp, 7.0_dp, &
      "CappedProtection", hedge, 0.5_dp)
    call assert_close("hedge amount conservation", &
      hedge%protected_amount + hedge%unprotected_amount, 14.0_dp, &
      1.0e-12_dp, failures)
    call assert_close("hedge RWA consistency", hedge%fx_rwa + hedge%cds_rwa, &
      hedge%unprotected_amount * 0.2_dp + &
      hedge%protected_amount * 0.3_dp, 1.0e-12_dp, failures)

    trades(1)%mtm = -100.0_dp
    agreement%values_type = "Percentage"
    agreement%threshold_counterparty = 0.10_dp
    agreement%minimum_transfer_counterparty = 0.05_dp
    agreement%initial_margin_counterparty = 0.02_dp
    call calc_rc(trades, rc, csa=agreement)
    call assert_close("percentage CSA RC", rc%rc, 13.0_dp, &
      1.0e-12_dp, failures)
  end subroutine test_modes_and_hedge

  subroutine test_sold_options(failures)
    integer, intent(inout) :: failures
    type(trade_t) :: trades(1)
    type(portfolio_result_t) :: result

    call trades(1)%configure_class("IRDSwaption")
    trades(1)%trade_type = "Option"
    trades(1)%buy_sell = "Sell"
    call calculate_portfolio(trades, result)
    call assert_integer("sold-option exposure count", size(result%exposures), &
      0, failures)
    call assert_close("sold-option EAD", result%total_ead, 0.0_dp, &
      0.0_dp, failures)
  end subroutine test_sold_options

  subroutine assert_close(name, actual, expected, tolerance, failures)
    character(len=*), intent(in) :: name
    real(dp), intent(in) :: actual
    real(dp), intent(in) :: expected
    real(dp), intent(in) :: tolerance
    integer, intent(inout) :: failures

    if (abs(actual - expected) > tolerance) then
      failures = failures + 1
      write(*, '(a,2(1x,es24.16))') "FAIL " // trim(name), actual, expected
    end if
  end subroutine assert_close

  subroutine assert_integer(name, actual, expected, failures)
    character(len=*), intent(in) :: name
    integer, intent(in) :: actual
    integer, intent(in) :: expected
    integer, intent(inout) :: failures

    if (actual /= expected) then
      failures = failures + 1
      write(*, '(a,2(1x,i0))') "FAIL " // trim(name), actual, expected
    end if
  end subroutine assert_integer

  subroutine assert_character(name, actual, expected, failures)
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: actual
    character(len=*), intent(in) :: expected
    integer, intent(inout) :: failures

    if (trim(actual) /= trim(expected)) then
      failures = failures + 1
      write(*, '(a,2(1x,a))') "FAIL " // trim(name), trim(actual), trim(expected)
    end if
  end subroutine assert_character

  subroutine assert_true(name, condition, failures)
    character(len=*), intent(in) :: name
    logical, intent(in) :: condition
    integer, intent(inout) :: failures

    if (.not. condition) then
      failures = failures + 1
      write(*, '(a)') "FAIL " // trim(name)
    end if
  end subroutine assert_true

  pure logical function starts_with(text, prefix) result(value)
    character(len=*), intent(in) :: text
    character(len=*), intent(in) :: prefix

    if (len_trim(text) < len_trim(prefix)) then
      value = .false.
    else
      value = text(1:len_trim(prefix)) == prefix(1:len_trim(prefix))
    end if
  end function starts_with

end program test_saccr
