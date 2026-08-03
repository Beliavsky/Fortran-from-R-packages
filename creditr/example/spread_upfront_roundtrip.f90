! SPDX-License-Identifier: GPL-3.0-only AND LicenseRef-ISDA-CDS-Standard-Model
program spread_upfront_roundtrip
  use creditr
  implicit none
  type(rate_quote_t), allocatable :: quotes(:)
  type(conventions_t) :: conventions
  type(zero_curve_t) :: curve
  type(cds_dates_t) :: dates
  integer :: status
  real(kind=dp) :: upfront, recovered_spread

  call read_rate_quotes_csv('data/usd_2014_04_15.csv', quotes, status)
  call add_conventions('USD', conventions, status)
  call build_zero_curve(make_date(2014, 4, 17), quotes, conventions, curve, status)
  call add_dates(make_date(2014, 4, 15), 'USD', dates, maturity=make_date(2019, 6, 20), status=status)

  upfront = spread_to_upfront(curve, dates, 1737.7289_dp, 500.0_dp, 0.4_dp, 1.0e7_dp)
  recovered_spread = upfront_to_spread(curve, dates, upfront, 500.0_dp, 0.4_dp, 1.0e7_dp)

  print '(a,f14.2)', 'Dirty upfront:   ', upfront
  print '(a,f14.6)', 'Recovered spread:', recovered_spread
end program spread_upfront_roundtrip
