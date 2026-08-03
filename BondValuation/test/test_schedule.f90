! SPDX-License-Identifier: GPL-3.0-only
program test_schedule
  use bondvaluation
  implicit none
  type(bond_terms) :: terms
  type(bond_schedule) :: schedule
  integer :: i

  terms%issue_date = date_from_ymd(2020, 1, 15)
  terms%maturity_date = date_from_ymd(2030, 1, 15)
  terms%coupon_frequency = 2
  terms%redemption_value = 100.0_dp
  terms%coupon_rate_percent = 5.0_dp
  terms%day_count_convention = dcc_act_act_icma
  call build_bond_schedule(terms, schedule)
  call assert_equal_int(schedule%status, 0, "regular schedule status")
  call assert_equal_int(size(schedule%coupon_dates), 20, "regular coupon count")
  do i = 1, size(schedule%coupon_payments)
    call assert_close(schedule%coupon_payments(i), 2.5_dp, 2.0e-13_dp, &
      "regular coupon")
  end do
  call assert_date(schedule%coupon_dates(1), 2020, 7, 15, "first regular coupon")
  call assert_date(schedule%coupon_dates(20), 2030, 1, 15, "maturity coupon")
  if (trim(schedule%first_coupon_type) /= "regular") error stop "regular first type"
  if (trim(schedule%last_coupon_type) /= "regular") error stop "regular last type"

  terms = bond_terms()
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
  call assert_equal_int(schedule%status, 0, "odd schedule status")
  if (.not. schedule%end_of_month_used) error stop "EOM inference failed"
  call assert_equal_int(size(schedule%coupon_dates), 12, "odd coupon count")
  call assert_date(schedule%coupon_dates(1), 2015, 2, 28, "odd first coupon")
  call assert_date(schedule%coupon_dates(11), 2020, 2, 29, "odd last regular coupon")
  call assert_date(schedule%coupon_dates(12), 2021, 4, 21, "odd maturity")
  if (trim(schedule%first_coupon_type) /= "long") error stop "odd first type"
  if (trim(schedule%last_coupon_type) /= "long") error stop "odd last type"
  call assert_close(schedule%coupon_payments(1), 6.555248618784531_dp, &
    2.0e-12_dp, "odd first coupon payment")
  do i = 2, 11
    call assert_close(schedule%coupon_payments(i), 2.625_dp, 2.0e-12_dp, &
      "odd regular coupon payment")
  end do
  call assert_close(schedule%coupon_payments(12), 5.991847826086961_dp, &
    2.0e-12_dp, "odd final coupon payment")
  print '(a)', "test_schedule: PASS"

contains

  subroutine assert_close(actual, expected, tolerance, label)
    real(dp), intent(in) :: actual, expected, tolerance
    character(len=*), intent(in) :: label
    if (abs(actual - expected) > tolerance) then
      write(*, '(a,2(1x,es24.16))') trim(label)//" mismatch:", actual, expected
      error stop 1
    end if
  end subroutine assert_close

  subroutine assert_equal_int(actual, expected, label)
    integer, intent(in) :: actual, expected
    character(len=*), intent(in) :: label
    if (actual /= expected) then
      write(*, '(a,2(1x,i0))') trim(label)//" mismatch:", actual, expected
      error stop 1
    end if
  end subroutine assert_equal_int

  subroutine assert_date(actual, year, month, day, label)
    type(date_type), intent(in) :: actual
    integer, intent(in) :: year, month, day
    character(len=*), intent(in) :: label
    if (actual%year /= year .or. actual%month /= month .or. actual%day /= day) then
      write(*, '(a,1x,a)') trim(label)//" mismatch:", date_to_string(actual)
      error stop 1
    end if
  end subroutine assert_date

end program test_schedule
