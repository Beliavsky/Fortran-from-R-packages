! SPDX-License-Identifier: GPL-3.0-only
module test_compat_callbacks
  use bondvaluation_kinds, only: rk => dp
  implicit none
contains
  function square_minus_two(x) result(value)
    real(rk), intent(in) :: x
    real(rk) :: value
    value = x * x - 2.0_rk
  end function square_minus_two

  function twice_x(x) result(value)
    real(rk), intent(in) :: x
    real(rk) :: value
    value = 2.0_rk * x
  end function twice_x
end module test_compat_callbacks

program test_compat
  use bondvaluation_kinds, only: rk => dp
  use bondvaluation_dates, only: date_type, date_from_ymd, date_to_string
  use bondvaluation_daycount, only: daycount_result
  use bondvaluation_schedule, only: bond_terms, bond_schedule
  use bondvaluation_pricing, only: bond_value_result, accrued_interest_result
  use bondvaluation_compat
  use test_compat_callbacks, only: square_minus_two, twice_x
  implicit none
  type(date_type) :: date, previous, following
  type(daycount_result) :: dc
  type(accrued_interest_result) :: ai
  type(bond_terms) :: terms
  type(bond_schedule) :: schedule
  type(bond_value_result) :: price_result, yield_result
  type(newton_result) :: root_result
  real(rk), allocatable :: payments(:)
  type(date_type) :: starts(2), ends(2)

  call assert_equal_int(leap(2024), 1, "leap")
  call assert_equal_int(DaysInMonth(date_from_ymd(2024, 2, 1)), 29, "days in month")
  call assert_equal_int(DaysInYear(2023), 365, "days in year")
  call assert_equal_int(DayDiff(date_from_ymd(2024, 1, 1), &
    date_from_ymd(2024, 1, 31)), 30, "day diff")
  date = Date_LDM(date_from_ymd(2024, 2, 10))
  call assert_date(date, 2024, 2, 29, "last day")
  call assert_equal_int(LDM(date), 1, "last-day indicator")
  call assert_equal_int(LeapDayInside(date_from_ymd(2024, 1, 1), &
    date_from_ymd(2024, 3, 1)), 1, "leap day inside")
  call assert_close(sumC([1.0_rk, 2.0_rk, 3.0_rk]), 6.0_rk, 0.0_rk, "sumC")
  call assert_equal_int(FirstMatch(4, [2, 4, 4]), 2, "FirstMatch")
  call assert_equal_int(FirstMatch(9, [2, 4, 4]), 4, "FirstMatch missing")

  dc = DIST(10, date_from_ymd(2024, 1, 1), date_from_ymd(2024, 7, 1))
  call assert_equal_int(dc%days_accrued, 182, "DIST days")
  call assert_close(dc%fraction, 182.0_rk / 365.0_rk, 1.0e-15_rk, "DIST fraction")
  ai = AccrInt(date_from_ymd(2024, 1, 1), date_from_ymd(2024, 7, 1), &
    5.0_rk, 10, 100.0_rk)
  call assert_close(ai%accrued_interest, 100.0_rk * 0.05_rk * 182.0_rk / 365.0_rk, &
    1.0e-14_rk, "AccrInt")

  starts = [date_from_ymd(2024, 1, 1), date_from_ymd(2024, 7, 1)]
  ends = [date_from_ymd(2024, 7, 1), date_from_ymd(2025, 1, 1)]
  payments = PayCalc(100.0_rk, 0.05_rk, 10, starts, ends)
  call assert_close(payments(1), 100.0_rk * 0.05_rk * 182.0_rk / 365.0_rk, &
    1.0e-14_rk, "PayCalc first")
  call assert_close(payments(2), 100.0_rk * 0.05_rk * 184.0_rk / 365.0_rk, &
    1.0e-14_rk, "PayCalc second")

  date = NumToDate(1, date_from_ymd(1970, 1, 1))
  call assert_date(date, 1970, 1, 2, "NumToDate")
  previous = CppPrevDate(date_from_ymd(2024, 8, 31), &
    date_from_ymd(2024, 1, 1), date_from_ymd(2025, 2, 28), 2, .true.)
  call assert_date(previous, 2024, 2, 29, "CppPrevDate")
  following = CppSuccDate(date_from_ymd(2024, 2, 29), &
    date_from_ymd(2025, 1, 1), date_from_ymd(2025, 2, 28), 2, .true.)
  call assert_date(following, 2024, 8, 31, "CppSuccDate")

  root_result = NewtonRaphson(square_minus_two, twice_x, 1.0_rk, 1.0e-13_rk)
  call assert_equal_int(root_result%status, 0, "Newton status")
  call assert_close(root_result%root, sqrt(2.0_rk), 2.0e-13_rk, "Newton root")
  call assert_close(dm_MyPriceEqn(1.025_rk, 1, 1.5_rk, 2.625_rk, &
    102.625_rk, 0.4_rk, 8.0_rk, 0.6_rk, 2), &
    -404.4253185365378_rk, 2.0e-11_rk, "price derivative")
  call assert_close(dm_MyPriceEqn(1.025_rk, 2, 1.5_rk, 2.625_rk, &
    102.625_rk, 0.4_rk, 8.0_rk, 0.6_rk, 2), &
    1905.349184218272_rk, 2.0e-10_rk, "price second derivative")
  call assert_close(ModDUR(1.025_rk, 1.5_rk, 2.625_rk, 102.625_rk, &
    0.4_rk, 8.0_rk, 0.6_rk, 2, 102.2966287255355_rk), &
    3.9534569572338527_rk, 2.0e-12_rk, "ModDUR")
  call assert_close(CONV(1.025_rk, 1.5_rk, 2.625_rk, 102.625_rk, &
    0.4_rk, 8.0_rk, 0.6_rk, 2, 102.2966287255355_rk), &
    9.312864011043446_rk, 2.0e-12_rk, "CONV")

  terms%issue_date = date_from_ymd(2020, 1, 15)
  terms%maturity_date = date_from_ymd(2030, 1, 15)
  terms%coupon_frequency = 2
  terms%redemption_value = 100.0_rk
  terms%coupon_rate_percent = 5.0_rk
  terms%day_count_convention = 2
  schedule = AnnivDates(terms)
  price_result = BondVal_Price(4.0_rk, date_from_ymd(2024, 4, 15), &
    terms, schedule=schedule)
  yield_result = BondVal_Yield(price_result%clean_price, &
    date_from_ymd(2024, 4, 15), terms, schedule=schedule)
  call assert_close(yield_result%yield_percent, 4.0_rk, 2.0e-10_rk, &
    "compat yield inversion")
  print '(a)', "test_compat: PASS"

contains

  subroutine assert_close(actual, expected, tolerance, label)
    real(rk), intent(in) :: actual, expected, tolerance
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

end program test_compat
