program xva_demo
  use xva
  implicit none

  type(trade_t) :: trades(2)
  type(csa_t) :: csa
  type(collateral_t) :: collateral
  type(simulation_data_t) :: sim_data
  type(regulatory_data_t) :: reg_data
  type(curve_t) :: credit_curve_counterparty
  type(curve_t) :: credit_curve_processing_organization
  type(curve_t) :: funding_curve
  type(curve_t) :: spot_rates
  type(supervisory_cva_data_t) :: supervisory
  type(xva_result_t) :: result

  call make_trades(trades)
  call make_curves(credit_curve_processing_organization, &
    credit_curve_counterparty, funding_curve)
  call load_curve_csv("data/spot_rates.csv", spot_rates)
  spot_rates%rates = spot_rates%rates / 10000.0_dp
  call load_supervisory_cva_data("data", supervisory)

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
  sim_data%num_simulations = 50
  sim_data%mean_reversion_a = 0.001_dp
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

  call xva_calculator(trades, csa, collateral, sim_data, reg_data, &
    credit_curve_processing_organization, credit_curve_counterparty, &
    funding_curve, spot_rates, 0.45_dp, 0.45_dp, .false., result, &
    seed=20250523, supervisory_data=supervisory)

  print '(a,es14.6)', "CVA simulated:       ", result%cva_simulated
  print '(a,es14.6)', "DVA simulated:       ", result%dva_simulated
  print '(a,es14.6)', "FCA simulated:       ", result%fca_simulated
  print '(a,es14.6)', "FBA simulated:       ", result%fba_simulated
  print '(a,es14.6)', "MVA simulated:       ", result%mva_simulated
  print '(a,es14.6)', "CVA SA-CCR:          ", result%cva_saccr
  print '(a,es14.6)', "DVA SA-CCR:          ", result%dva_saccr
  print '(a,es14.6)', "SA-CVA IR charge:    ", &
    result%cva_capital%interest_rate_charge
  print '(a,es14.6)', "SA-CVA CS charge:    ", &
    result%cva_capital%credit_spread_charge
  print '(a,es14.6)', "CVA capital total:   ", result%cva_capital%total_charge
  print '(a,es14.6)', "KVA:                 ", result%kva
  print '(a,es14.6)', "Regulatory EAD:      ", result%ead
  print '(a,es14.6)', "Effective maturity:  ", result%effective_maturity

contains

  subroutine make_trades(values)
    type(trade_t), intent(out) :: values(2)

    call values(1)%configure_class("IRDSwap")
    values(1)%external_id = "ext1"
    values(1)%notional = 1.0_dp
    values(1)%mtm = 0.045_dp
    values(1)%currency = "USD"
    values(1)%si = 0.0_dp
    values(1)%ei = 7.0_dp
    values(1)%buy_sell = "Sell"
    values(1)%pay_leg_rate = 0.05_dp

    call values(2)%configure_class("IRDSwap")
    values(2)%external_id = "ext2"
    values(2)%notional = 1.0_dp
    values(2)%mtm = -0.065_dp
    values(2)%currency = "USD"
    values(2)%si = 0.0_dp
    values(2)%ei = 10.0_dp
    values(2)%buy_sell = "Buy"
    values(2)%pay_leg_rate = 0.05_dp
  end subroutine make_trades

  subroutine make_curves(processing_organization, counterparty, funding)
    type(curve_t), intent(out) :: processing_organization
    type(curve_t), intent(out) :: counterparty
    type(curve_t), intent(out) :: funding
    real(dp), parameter :: tenors(7) = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, &
      5.0_dp, 6.0_dp, 10.0_dp]

    allocate(counterparty%tenors(7), counterparty%rates(7), &
      processing_organization%tenors(7), processing_organization%rates(7), &
      funding%tenors(7), funding%rates(7))
    counterparty%tenors = tenors
    counterparty%rates = [3.0_dp, 10.0_dp, 20.0_dp, 40.0_dp, 66.0_dp, &
      99.0_dp, 150.0_dp]
    processing_organization%tenors = tenors
    processing_organization%rates = [4.0_dp, 11.0_dp, 23.0_dp, 47.0_dp, &
      76.0_dp, 110.0_dp, 160.0_dp]
    funding%tenors = tenors
    funding%rates = [4.0_dp, 17.0_dp, 43.0_dp, 47.0_dp, 76.0_dp, &
      90.0_dp, 110.0_dp]
  end subroutine make_curves

end program xva_demo
