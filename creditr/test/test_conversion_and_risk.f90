! SPDX-License-Identifier: GPL-3.0-only AND LicenseRef-ISDA-CDS-Standard-Model
program test_conversion_and_risk
  use creditr
  implicit none
  type(rate_quote_t), allocatable :: quotes(:)
  type(conventions_t) :: conventions
  type(zero_curve_t) :: curve
  type(cds_dates_t) :: dates
  real(kind=dp) :: upfront, spread, dv01, rr01, ten_percent
  integer :: status, failures

  failures = 0
  call read_rate_quotes_csv('data/usd_2014_04_15.csv', quotes, status)
  call add_conventions('USD', conventions, status)
  call build_zero_curve(make_date(2014, 4, 17), quotes, conventions, curve, status)
  call add_dates(make_date(2014, 4, 15), 'USD', dates, maturity=make_date(2019, 6, 20), status=status)

  upfront = spread_to_upfront(curve, dates, 243.28_dp, 100.0_dp, 0.4_dp, 1.0e7_dp, .false., status)
  spread = upfront_to_spread(curve, dates, upfront, 100.0_dp, 0.4_dp, 1.0e7_dp, .false., status)
  call check_close(spread, 243.28_dp, 2.0e-7_dp, 'upfront/spread inversion')

  dv01 = spread_dv01(curve, dates, 243.28_dp, 100.0_dp, 0.4_dp, 1.0e7_dp)
  rr01 = rec_risk_01(curve, dates, 243.28_dp, 100.0_dp, 0.4_dp, 1.0e7_dp)
  ten_percent = cs10(curve, dates, 243.28_dp, 100.0_dp, 0.4_dp, 1.0e7_dp)
  call check_close(dv01, 4282.0492767_dp, 2.0e-3_dp, 'spread DV01')
  call check_close(rr01, -1109.6460261_dp, 2.0e-3_dp, 'recovery risk 01')
  call check(ten_percent > 20.0_dp * dv01, 'CS10 exceeds 20 spread DV01 increments')

  if (failures /= 0) error stop 'test_conversion_and_risk failed'
  print '(a)', 'test_conversion_and_risk: PASS'

contains

  subroutine check(condition, label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    if (.not. condition) then
      failures = failures + 1
      print '(a)', 'FAIL: ' // trim(label)
    end if
  end subroutine check

  subroutine check_close(actual, expected, tolerance, label)
    real(kind=dp), intent(in) :: actual, expected, tolerance
    character(len=*), intent(in) :: label
    call check(abs(actual - expected) <= tolerance, label)
  end subroutine check_close

end program test_conversion_and_risk
