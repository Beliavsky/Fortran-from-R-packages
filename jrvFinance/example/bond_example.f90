! SPDX-License-Identifier: GPL-2.0-or-later
program bond_example
  use jrvfinance, only: dp, date_t, parse_date, bond_price, bond_yield, &
    bond_duration, coupons_dates, date_string
  implicit none
  type(date_t) :: settle, mature
  type(date_t), allocatable :: dates(:)
  real(dp) :: price, yield_rate

  settle = parse_date('2012-04-15')
  mature = parse_date('2022-01-01')
  price = bond_price(settle, mature, 0.08_dp, 0.088843_dp)
  yield_rate = bond_yield(settle, mature, 0.08_dp, price)
  dates = coupons_dates(settle, mature)
  print '(a,f10.4)', 'Clean price: ', price
  print '(a,f10.6)', 'Yield: ', yield_rate
  print '(a,f10.4)', 'Macaulay duration: ', &
    bond_duration(settle, mature, 0.08_dp, yield_rate)
  print '(a,i0)', 'Coupons remaining: ', size(dates)
  print '(a,a)', 'Next coupon: ', date_string(dates(1))
end program bond_example
