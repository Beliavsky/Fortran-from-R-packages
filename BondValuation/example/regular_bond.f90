! SPDX-License-Identifier: GPL-3.0-only
program regular_bond
  use bondvaluation
  implicit none
  type(bond_terms) :: terms
  type(bond_value_result) :: value, recovered

  terms%issue_date = date_from_ymd(2020, 1, 15)
  terms%maturity_date = date_from_ymd(2030, 1, 15)
  terms%coupon_frequency = 2
  terms%redemption_value = 100.0_dp
  terms%coupon_rate_percent = 5.0_dp
  terms%day_count_convention = dcc_act_act_icma

  value = bond_price(4.0_dp, date_from_ymd(2024, 4, 15), terms)
  recovered = bond_yield(value%clean_price, date_from_ymd(2024, 4, 15), terms)
  print '(a,f12.6)', "Clean price:  ", value%clean_price
  print '(a,f12.6)', "Dirty price:  ", value%dirty_price
  print '(a,f12.6)', "Accrued:      ", value%accrued_interest
  print '(a,f12.6)', "Recovered YTM:", recovered%yield_percent
end program regular_bond
