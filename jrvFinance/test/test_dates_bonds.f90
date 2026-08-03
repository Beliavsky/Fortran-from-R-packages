! SPDX-License-Identifier: GPL-2.0-or-later
program test_dates_bonds
  use jrvfinance, only: dp, date_t, date, parse_date, date_string, edate, &
    daycount_actual, daycount_30_360, year_fraction, coupons_n, &
    coupons_next, coupons_prev, coupons_dates, bond_cashflows, bond_tcf, &
    bond_price, bond_yield, bond_duration, JRV_OK
  implicit none
  type(date_t) :: settle, mature, d
  type(date_t), allocatable :: dates(:)
  type(bond_cashflows) :: tcf
  real(dp) :: price, ytm, dur
  integer :: status

  settle = parse_date('2012-04-15', status)
  mature = parse_date('2022-01-01', status)
  call check(status == JRV_OK, 'parse date')
  d = edate(parse_date('2007-02-28'), 4)
  call check(date_string(d) == '2007-06-30', 'edate end of month')
  d = edate(parse_date('2005-05-17'), -8)
  call check(date_string(d) == '2004-09-17', 'edate negative')
  call check(daycount_actual(date(2020,2,28),date(2020,3,1)) == 2, 'actual day count')
  call check(daycount_30_360(date(2020,1,30),date(2020,2,28)) == 28, '30/360 day count')
  call check(abs(year_fraction(date(2020,1,1),date(2020,7,1), &
    date(2020,1,1),date(2020,7,1),2,'ACT/ACT')-0.5_dp) < 1.0e-12_dp, &
    'actual actual fraction')

  call check(coupons_n(settle,mature,2) == 20, 'coupon count')
  call check(date_string(coupons_next(settle,mature,2)) == '2012-07-01', 'next coupon')
  call check(date_string(coupons_prev(settle,mature,2)) == '2012-01-01', 'previous coupon')
  dates = coupons_dates(settle,mature,2)
  call check(size(dates) == 20 .and. date_string(dates(size(dates))) == &
    '2022-01-01', 'coupon dates')

  tcf = bond_tcf(settle,mature,0.08_dp,2,'30/360',100.0_dp)
  call check(tcf%status == JRV_OK .and. size(tcf%time) == 20, 'bond cash flows')
  call check(abs(tcf%cashflow(size(tcf%cashflow))-104.0_dp) < 1.0e-12_dp, &
    'redemption cash flow')

  price = bond_price(settle,mature,0.08_dp,0.088843_dp,2,'30/360',2.0_dp,100.0_dp,status)
  call check(status == JRV_OK .and. price > 90.0_dp .and. price < 100.0_dp, 'bond price')
  ytm = bond_yield(settle,mature,0.08_dp,price,2,'30/360',2.0_dp,100.0_dp,status)
  call check(status == JRV_OK .and. abs(ytm-0.088843_dp) < 1.0e-7_dp, 'bond yield round trip')
  dur = bond_duration(settle,mature,0.08_dp,ytm,2,'30/360',.false.,2.0_dp,100.0_dp,status)
  call check(status == JRV_OK .and. dur > 5.0_dp .and. dur < 8.0_dp, 'bond duration')

  print '(a)', 'test_dates_bonds: PASS'
contains
  subroutine check(condition, label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    if (.not. condition) then
      write(*,'(a)') 'FAIL: '//trim(label)
      error stop 1
    end if
  end subroutine check
end program test_dates_bonds
