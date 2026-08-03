! SPDX-License-Identifier: GPL-3.0-only
program day_count_comparison
  use bondvaluation
  implicit none
  type(accrued_interest_result) :: result
  integer :: dcc

  print '(a)', "DCC  convention                           days       fraction       accrued"
  do dcc = 1, 16
    result = accrued_interest(date_from_ymd(2011, 8, 31), &
      date_from_ymd(2012, 2, 29), 5.25_dp, dcc, 10000.0_dp, 2, &
      date_from_ymd(2021, 8, 31), 2012, .true.)
    print '(i3,2x,a36,2x,i4,2x,f13.9,2x,f12.6)', dcc, day_count_name(dcc), &
      result%days_accrued, result%year_fraction_value, result%accrued_interest
  end do
end program day_count_comparison
