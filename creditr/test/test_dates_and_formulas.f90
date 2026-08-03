! SPDX-License-Identifier: GPL-3.0-only AND LicenseRef-ISDA-CDS-Standard-Model
program test_dates_and_formulas
  use creditr
  implicit none
  type(cds_dates_t) :: dates
  type(conventions_t) :: conventions
  integer :: status, failures
  real(kind=dp) :: pd, spread

  failures = 0
  call add_dates(make_date(2014, 4, 15), 'USD', dates, maturity=make_date(2019, 6, 20), status=status)
  call check(status == creditr_ok, 'add_dates status')
  call check_date(dates%stepin_date, 2014, 4, 16, 'step-in date')
  call check_date(dates%value_date, 2014, 4, 18, 'value date')
  call check_date(dates%start_date, 2014, 3, 20, 'start date')
  call check_date(dates%first_coupon_date, 2014, 6, 20, 'first coupon')
  call check_date(dates%penultimate_coupon_date, 2019, 3, 20, 'penultimate coupon')
  call check_date(dates%base_date, 2014, 4, 17, 'base date')
  call check_date(dates%backstop_date, 2014, 2, 14, 'backstop date')

  call add_dates(make_date(2015, 9, 20), 'JPY', dates, tenor_years=5, status=status)
  call check_date(dates%base_date, 2015, 9, 25, 'JPY holiday base date')

  call add_dates(make_date(2010, 6, 18), 'USD', dates, tenor_years=5, status=status)
  call check_date(dates%end_date, 2015, 6, 20, 'weekend maturity remains unadjusted')

  call add_conventions('USD', conventions, status)
  call check(conventions%fixed_frequency_months == 6, 'USD fixed frequency')
  call check(trim(conventions%fixed_dcc) == '30/360', 'USD fixed DCC')
  call add_conventions('EUR', conventions, status)
  call check(conventions%fixed_frequency_months == 12, 'EUR fixed frequency')
  call add_conventions('JPY', conventions, status)
  call check(trim(conventions%mm_dcc) == 'ACT/365', 'JPY money-market DCC')

  pd = spread_to_pd(250.0_dp, 0.4_dp, 5.0_dp)
  spread = pd_to_spread(pd, 0.4_dp, 5.0_dp)
  call check_close(spread, 250.0_dp, 1.0e-11_dp, 'spread/PD inversion')
  call check_close(implied_rr(0.1_dp, 200.0_dp, 5.0_dp), 5.087784189700955_dp, 1.0e-8_dp, 'implied recovery')
  call check_close(pv01(500000.0_dp, 1.0e7_dp, 500.0_dp, 100.0_dp), 1.25_dp, 1.0e-13_dp, 'PV01')

  if (failures /= 0) error stop 'test_dates_and_formulas failed'
  print '(a)', 'test_dates_and_formulas: PASS'

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

  subroutine check_date(actual, year, month, day, label)
    type(date_t), intent(in) :: actual
    integer, intent(in) :: year, month, day
    character(len=*), intent(in) :: label
    call check(actual%year == year .and. actual%month == month .and. actual%day == day, label)
  end subroutine check_date

end program test_dates_and_formulas
