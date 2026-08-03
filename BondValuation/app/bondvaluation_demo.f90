! SPDX-License-Identifier: GPL-3.0-only
program bondvaluation_demo
  use bondvaluation
  implicit none
  type(bond_terms) :: terms
  type(bond_schedule) :: schedule
  type(bond_value_result) :: value, recovered
  integer :: i

  terms%issue_date = date_from_ymd(2013, 11, 30)
  terms%maturity_date = date_from_ymd(2021, 4, 21)
  terms%coupon_frequency = 2
  terms%first_coupon_date = date_from_ymd(2015, 2, 28)
  terms%last_coupon_date = date_from_ymd(2020, 2, 29)
  terms%first_interest_accrual_date = terms%issue_date
  terms%redemption_value = 100.0_dp
  terms%coupon_rate_percent = 5.25_dp
  terms%day_count_convention = dcc_act_act_icma

  call build_bond_schedule(terms, schedule)
  if (schedule%status /= 0) error stop "could not build schedule"
  print '(a)', "Coupon schedule"
  print '(a)', "date        coupon"
  do i = 1, size(schedule%coupon_dates)
    print '(a,2x,f10.6)', date_to_string(schedule%coupon_dates(i)), &
      schedule%coupon_payments(i)
  end do

  value = bond_price(5.0_dp, date_from_ymd(2014, 10, 15), terms, &
    supplied_schedule=schedule)
  recovered = bond_yield(value%clean_price, date_from_ymd(2014, 10, 15), &
    terms, supplied_schedule=schedule)
  print '(/,a,f12.6)', "Clean price:       ", value%clean_price
  print '(a,f12.6)',   "Accrued interest:  ", value%accrued_interest
  print '(a,f12.6)',   "Dirty price:       ", value%dirty_price
  print '(a,f12.6)',   "Recovered YTM (%): ", recovered%yield_percent
  print '(a,f12.6)',   "Modified duration: ", value%modified_duration_years
  print '(a,f12.6)',   "Convexity:         ", value%convexity_years
end program bondvaluation_demo
