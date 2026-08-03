! SPDX-License-Identifier: GPL-3.0-only AND LicenseRef-ISDA-CDS-Standard-Model
program creditr_demo
  use creditr
  implicit none
  type(rate_quote_t), allocatable :: quotes(:)
  type(conventions_t) :: conventions
  type(zero_curve_t) :: curve
  type(cds_contract_t) :: contract
  type(cds_result_t) :: result
  integer :: status

  call read_rate_quotes_csv('data/usd_2014_04_15.csv', quotes, status)
  call add_conventions('USD', conventions, status)
  call build_zero_curve(make_date(2014, 4, 17), quotes, conventions, curve, status)

  contract%name = 'Example reference entity'
  contract%trade_date = make_date(2014, 4, 15)
  contract%maturity = make_date(2019, 6, 20)
  contract%use_maturity = .true.
  contract%spread_bps = 243.28_dp
  contract%coupon_bps = 100.0_dp
  contract%recovery = 0.4_dp
  contract%currency = 'USD'
  contract%notional = 1.0e7_dp

  result = cds(contract, curve, status, quotes)
  if (status /= creditr_ok) error stop 'CDS valuation failed'

  print '(a)', 'Based on the ISDA CDS Standard Model (version 1.0),'
  print '(a)', 'developed and supported in collaboration with Markit Group Ltd.'
  print '(a)'
  print '(a,f14.2)', 'Clean principal:    ', result%principal
  print '(a,f14.2)', 'Dirty upfront:      ', result%upfront
  print '(a,f14.6)', 'Clean price:        ', result%price
  print '(a,f14.6)', 'Hazard rate:        ', result%hazard_rate
  print '(a,f14.2)', 'Spread DV01:        ', result%spread_dv01
  print '(a,f14.2)', 'Interest-rate DV01: ', result%ir_dv01
  print '(a,f14.2)', 'Recovery risk 01:   ', result%recovery_risk_01
end program creditr_demo
