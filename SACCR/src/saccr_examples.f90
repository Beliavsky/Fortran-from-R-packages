module saccr_examples
  use trading, only : dp, trade_t, csa_t, collateral_t
  use saccr_types, only : portfolio_result_t, exposure_result_t, fx_hedge_result_t
  use saccr_portfolio, only : calculate_portfolio
  use saccr_core, only : apply_fx_hedge
  implicit none
  private

  public :: example_ird
  public :: example_fx
  public :: example_credit
  public :: example_commodity
  public :: example_ird_credit
  public :: example_ird_commodity_margined
  public :: example_basis_volatility
  public :: example_fx_hedge

contains

  subroutine example_ird(result)
    type(portfolio_result_t), intent(out) :: result
    type(trade_t) :: trades(3)

    call trades(1)%configure_class("IRDSwap")
    trades(1)%external_id = "ext_1"
    trades(1)%notional = 10000.0_dp
    trades(1)%mtm = 30.0_dp
    trades(1)%currency = "USD"
    trades(1)%si = 0.0_dp
    trades(1)%ei = 10.0_dp
    trades(1)%buy_sell = "Buy"
    trades(1)%counterparty = "IRD Example"

    call trades(2)%configure_class("IRDSwap")
    trades(2)%external_id = "ext_2"
    trades(2)%notional = 10000.0_dp
    trades(2)%mtm = -20.0_dp
    trades(2)%currency = "USD"
    trades(2)%si = 0.0_dp
    trades(2)%ei = 4.0_dp
    trades(2)%buy_sell = "Sell"
    trades(2)%counterparty = "IRD Example"

    call trades(3)%configure_class("IRDSwaption")
    trades(3)%external_id = "ext_3"
    trades(3)%notional = 5000.0_dp
    trades(3)%mtm = 50.0_dp
    trades(3)%currency = "EUR"
    trades(3)%si = 1.0_dp
    trades(3)%ei = 11.0_dp
    trades(3)%buy_sell = "Buy"
    trades(3)%option_type = "Put"
    trades(3)%underlying_price = 0.06_dp
    trades(3)%strike_price = 0.05_dp
    trades(3)%counterparty = "IRD Example"

    call calculate_portfolio(trades, result)
  end subroutine example_ird

  subroutine example_fx(result)
    type(portfolio_result_t), intent(out) :: result
    type(trade_t) :: trades(3)

    call trades(1)%configure_class("FxForward")
    trades(1)%external_id = "ext_1"
    trades(1)%notional = 10000.0_dp
    trades(1)%mtm = 30.0_dp
    trades(1)%ccy_pair = "EUR/USD"
    trades(1)%si = 0.0_dp
    trades(1)%ei = 10.0_dp
    trades(1)%buy_sell = "Buy"

    call trades(2)%configure_class("FxForward")
    trades(2)%external_id = "ext_2"
    trades(2)%notional = 20000.0_dp
    trades(2)%mtm = -20.0_dp
    trades(2)%ccy_pair = "EUR/USD"
    trades(2)%si = 0.0_dp
    trades(2)%ei = 4.0_dp
    trades(2)%buy_sell = "Sell"

    call trades(3)%configure_class("FxForward")
    trades(3)%external_id = "ext_3"
    trades(3)%notional = 5000.0_dp
    trades(3)%mtm = 50.0_dp
    trades(3)%ccy_pair = "GBP/USD"
    trades(3)%si = 1.0_dp
    trades(3)%ei = 11.0_dp
    trades(3)%buy_sell = "Sell"

    call calculate_portfolio(trades, result)
  end subroutine example_fx

  subroutine example_credit(result)
    type(portfolio_result_t), intent(out) :: result
    type(trade_t) :: trades(3)

    call make_credit_examples(trades, "Credit Example")
    call calculate_portfolio(trades, result)
  end subroutine example_credit

  subroutine example_commodity(result)
    type(portfolio_result_t), intent(out) :: result
    type(trade_t) :: trades(3)

    call make_commodity_examples(trades, "Commodity Example", .false.)
    call calculate_portfolio(trades, result)
  end subroutine example_commodity

  subroutine example_ird_credit(result)
    type(portfolio_result_t), intent(out) :: result
    type(trade_t) :: credit(3)
    type(trade_t) :: ird(3)
    type(trade_t) :: trades(6)

    call make_credit_examples(credit, "IRDCredit Example")
    call make_ird_examples(ird, "IRDCredit Example")
    trades(:3) = credit
    trades(4:) = ird
    call calculate_portfolio(trades, result)
  end subroutine example_ird_credit

  subroutine example_ird_commodity_margined(result)
    type(portfolio_result_t), intent(out) :: result
    type(trade_t) :: commodity(3)
    type(trade_t) :: ird(3)
    type(trade_t) :: trades(6)
    type(csa_t) :: agreements(1)
    type(collateral_t) :: collateral(2)

    call make_commodity_examples(commodity, "IRDCommMargined", .true.)
    call make_ird_examples(ird, "IRDCommMargined")
    trades(:3) = commodity
    trades(4:) = ird

    agreements(1)%id = "csa_1"
    agreements(1)%threshold_counterparty = 0.0_dp
    agreements(1)%threshold_processing_organization = 0.0_dp
    agreements(1)%initial_margin_counterparty = 0.0_dp
    agreements(1)%initial_margin_processing_organization = 0.0_dp
    agreements(1)%minimum_transfer_counterparty = 5.0_dp
    agreements(1)%minimum_transfer_processing_organization = 5.0_dp
    agreements(1)%mpor_days = 10.0_dp
    agreements(1)%remargin_frequency = 5.0_dp
    agreements(1)%counterparty = "IRDCommMargined"
    agreements(1)%values_type = "Actual"
    allocate(agreements(1)%currencies(2), agreements(1)%trade_groups(2))
    agreements(1)%currencies = ["USD", "EUR"]
    agreements(1)%trade_groups = ["IRD      ", "Commodity"]

    collateral(1)%id = "coll_1"
    collateral(1)%amount = 150.0_dp
    collateral(1)%csa_id = "csa_1"
    collateral(1)%collateral_type = "ICA"
    collateral(2)%id = "coll_2"
    collateral(2)%amount = 50.0_dp
    collateral(2)%csa_id = "csa_1"
    collateral(2)%collateral_type = "VariationMargin"

    call calculate_portfolio(trades, result, csas=agreements, collaterals=collateral)
  end subroutine example_ird_commodity_margined

  subroutine example_basis_volatility(result)
    type(portfolio_result_t), intent(out) :: result
    type(trade_t) :: trades(4)

    call trades(1)%configure_class("IRDSwap")
    trades(1)%external_id = "ext_1"
    trades(1)%notional = 10000.0_dp
    trades(1)%mtm = 30.0_dp
    trades(1)%currency = "USD"
    trades(1)%si = 0.0_dp
    trades(1)%ei = 10.0_dp
    trades(1)%buy_sell = "Buy"
    trades(1)%pay_leg_type = "Float"
    trades(1)%pay_leg_ref = "CDOR"
    trades(1)%pay_leg_tenor = "1M"
    trades(1)%rec_leg_type = "Float"
    trades(1)%rec_leg_ref = "CORRA"
    trades(1)%rec_leg_tenor = "3M"

    call trades(2)%configure_class("CommSwap")
    trades(2)%external_id = "ext_2"
    trades(2)%notional = 10000.0_dp
    trades(2)%mtm = -20.0_dp
    trades(2)%currency = "USD"
    trades(2)%si = 0.0_dp
    trades(2)%ei = 4.0_dp
    trades(2)%buy_sell = "Sell"
    trades(2)%subclass = "Energy"
    trades(2)%commodity_type = "Oil/Gas"
    trades(2)%pay_leg_type = "Commodity"
    trades(2)%pay_leg_ref = "Brent"
    trades(2)%rec_leg_type = "Commodity"
    trades(2)%rec_leg_ref = "Gas"

    call trades(3)%configure_class("IRDSwapVol")
    trades(3)%external_id = "ext_3"
    trades(3)%notional = 5000.0_dp
    trades(3)%mtm = 50.0_dp
    trades(3)%currency = "EUR"
    trades(3)%si = 1.0_dp
    trades(3)%ei = 11.0_dp
    trades(3)%buy_sell = "Sell"
    trades(3)%underlying_instrument = "CDOR"
    trades(3)%vol_strike = 0.2_dp
    trades(3)%annualization_factor = 252.0_dp

    call trades(4)%configure_class("IRDSwap")
    trades(4)%external_id = "ext_4"
    trades(4)%notional = 10000.0_dp
    trades(4)%mtm = 30.0_dp
    trades(4)%currency = "USD"
    trades(4)%si = 0.0_dp
    trades(4)%ei = 10.0_dp
    trades(4)%buy_sell = "Buy"

    call calculate_portfolio(trades, result)
  end subroutine example_basis_volatility

  subroutine example_fx_hedge(exposure, hedge, rwa_fx_counterparty, &
      rwa_cds_counterparty, ead_cds, approach, protection_percentage)
    type(exposure_result_t), intent(out) :: exposure
    type(fx_hedge_result_t), intent(out) :: hedge
    real(dp), intent(in) :: rwa_fx_counterparty
    real(dp), intent(in) :: rwa_cds_counterparty
    real(dp), intent(in) :: ead_cds
    character(len=*), intent(in) :: approach
    real(dp), intent(in), optional :: protection_percentage
    type(portfolio_result_t) :: portfolio
    type(trade_t) :: trades(1)

    call trades(1)%configure_class("FxForward")
    trades(1)%external_id = "ext_1"
    trades(1)%notional = 250.0_dp
    trades(1)%mtm = 0.0_dp
    trades(1)%ccy_pair = "EUR/USD"
    trades(1)%si = 0.0_dp
    trades(1)%ei = 10.0_dp
    trades(1)%buy_sell = "Buy"

    call calculate_portfolio(trades, portfolio)
    exposure = portfolio%exposures(1)
    if (present(protection_percentage)) then
      call apply_fx_hedge(exposure%ead, exposure%addon%addon, &
        rwa_fx_counterparty, rwa_cds_counterparty, ead_cds, approach, hedge, &
        protection_percentage)
    else
      call apply_fx_hedge(exposure%ead, exposure%addon%addon, &
        rwa_fx_counterparty, rwa_cds_counterparty, ead_cds, approach, hedge)
    end if
  end subroutine example_fx_hedge

  subroutine make_ird_examples(trades, counterparty)
    type(trade_t), intent(out) :: trades(3)
    character(len=*), intent(in) :: counterparty

    call trades(1)%configure_class("IRDSwap")
    trades(1)%external_id = "ext_4"
    trades(1)%notional = 10000.0_dp
    trades(1)%mtm = 30.0_dp
    trades(1)%currency = "USD"
    trades(1)%si = 0.0_dp
    trades(1)%ei = 10.0_dp
    trades(1)%buy_sell = "Buy"
    trades(1)%counterparty = counterparty

    call trades(2)%configure_class("IRDSwap")
    trades(2)%external_id = "ext_5"
    trades(2)%notional = 10000.0_dp
    trades(2)%mtm = -20.0_dp
    trades(2)%currency = "USD"
    trades(2)%si = 0.0_dp
    trades(2)%ei = 4.0_dp
    trades(2)%buy_sell = "Sell"
    trades(2)%counterparty = counterparty

    call trades(3)%configure_class("IRDSwaption")
    trades(3)%external_id = "ext_6"
    trades(3)%notional = 5000.0_dp
    trades(3)%mtm = 50.0_dp
    trades(3)%currency = "EUR"
    trades(3)%si = 1.0_dp
    trades(3)%ei = 11.0_dp
    trades(3)%buy_sell = "Buy"
    trades(3)%option_type = "Put"
    trades(3)%underlying_price = 0.06_dp
    trades(3)%strike_price = 0.05_dp
    trades(3)%counterparty = counterparty
  end subroutine make_ird_examples

  subroutine make_credit_examples(trades, counterparty)
    type(trade_t), intent(out) :: trades(3)
    character(len=*), intent(in) :: counterparty

    call trades(1)%configure_class("CDS")
    trades(1)%external_id = "ext_1"
    trades(1)%notional = 10000.0_dp
    trades(1)%mtm = 20.0_dp
    trades(1)%currency = "USD"
    trades(1)%si = 0.0_dp
    trades(1)%ei = 3.0_dp
    trades(1)%buy_sell = "Buy"
    trades(1)%subclass = "AA"
    trades(1)%reference_entity = "FirmA"
    trades(1)%counterparty = counterparty

    call trades(2)%configure_class("CDS")
    trades(2)%external_id = "ext_2"
    trades(2)%notional = 10000.0_dp
    trades(2)%mtm = -40.0_dp
    trades(2)%currency = "EUR"
    trades(2)%si = 0.0_dp
    trades(2)%ei = 6.0_dp
    trades(2)%buy_sell = "Sell"
    trades(2)%subclass = "BBB"
    trades(2)%reference_entity = "FirmB"
    trades(2)%counterparty = counterparty

    call trades(3)%configure_class("CDX")
    trades(3)%external_id = "ext_3"
    trades(3)%notional = 10000.0_dp
    trades(3)%mtm = 0.0_dp
    trades(3)%currency = "USD"
    trades(3)%si = 0.0_dp
    trades(3)%ei = 5.0_dp
    trades(3)%buy_sell = "Buy"
    trades(3)%subclass = "IG"
    trades(3)%reference_entity = "CDX.IG"
    trades(3)%counterparty = counterparty
  end subroutine make_credit_examples

  subroutine make_commodity_examples(trades, counterparty, set_currency)
    type(trade_t), intent(out) :: trades(3)
    character(len=*), intent(in) :: counterparty
    logical, intent(in) :: set_currency

    call trades(1)%configure_class("Commodity")
    trades(1)%external_id = "ext_1"
    trades(1)%notional = 10000.0_dp
    trades(1)%mtm = -50.0_dp
    trades(1)%si = 0.0_dp
    trades(1)%ei = 0.75_dp
    trades(1)%buy_sell = "Buy"
    trades(1)%subclass = "Energy"
    trades(1)%commodity_type = "Oil/Gas"
    trades(1)%counterparty = counterparty

    call trades(2)%configure_class("Commodity")
    trades(2)%external_id = "ext_2"
    trades(2)%notional = 20000.0_dp
    trades(2)%mtm = -30.0_dp
    trades(2)%si = 0.0_dp
    trades(2)%ei = 2.0_dp
    trades(2)%buy_sell = "Sell"
    trades(2)%subclass = "Energy"
    trades(2)%commodity_type = "Oil/Gas"
    trades(2)%counterparty = counterparty

    call trades(3)%configure_class("Commodity")
    trades(3)%external_id = "ext_3"
    trades(3)%notional = 10000.0_dp
    trades(3)%mtm = 100.0_dp
    trades(3)%si = 0.0_dp
    trades(3)%ei = 5.0_dp
    trades(3)%buy_sell = "Buy"
    trades(3)%subclass = "Metals"
    trades(3)%commodity_type = "Silver"
    trades(3)%counterparty = counterparty

    if (set_currency) trades%currency = "USD"
  end subroutine make_commodity_examples

end module saccr_examples
