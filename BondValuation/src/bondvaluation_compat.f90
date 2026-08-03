! SPDX-License-Identifier: GPL-3.0-only
module bondvaluation_compat
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use bondvaluation_kinds, only: real_kind => dp
  use bondvaluation_dates, only: date_type, date_from_ymd, valid_date, &
    leap_year, days_in_month, days_in_year, day_diff, last_day_of_month, &
    is_last_day_of_month, date_to_serial, serial_to_date
  use bondvaluation_daycount, only: daycount_result, day_count_fraction, &
    leap_day_inside
  use bondvaluation_schedule, only: bond_terms, bond_schedule, build_bond_schedule
  use bondvaluation_pricing, only: accrued_interest_result, dirty_price_result, &
    bond_value_result, accrued_interest, dirty_price, bond_price, bond_yield
  implicit none
  private

  type, public :: newton_result
    real(real_kind) :: root = 0.0_real_kind
    integer :: iterations = 0
    real(real_kind) :: function_value = 0.0_real_kind
    integer :: status = 0
  end type newton_result

  abstract interface
    function scalar_function(x) result(value)
      import real_kind
      real(real_kind), intent(in) :: x
      real(real_kind) :: value
    end function scalar_function
  end interface

  public :: AccrInt, AnnivDates, DP, BondVal_Price, BondVal_Yield
  public :: leap, LDM, DaysInMonth, DaysInYear, DayDiff, Date_LDM
  public :: sumC, FirstMatch, LeapDayInside, DIST, PayCalc
  public :: NumToDate, CppPrevDate, CppSuccDate, NewtonRaphson
  public :: dm_MyPriceEqn, ModDUR, CONV

contains

  function AccrInt(start_date, end_date, coupon_rate_percent, dcc, &
                   redemption_value, coupon_frequency, maturity, &
                   next_coupon_year, end_of_month, anniversary_dates) result(out)
    type(date_type), intent(in) :: start_date, end_date
    real(real_kind), intent(in) :: coupon_rate_percent
    integer, intent(in) :: dcc
    real(real_kind), intent(in), optional :: redemption_value
    integer, intent(in), optional :: coupon_frequency
    type(date_type), intent(in), optional :: maturity
    integer, intent(in), optional :: next_coupon_year
    logical, intent(in), optional :: end_of_month
    type(date_type), intent(in), optional :: anniversary_dates(:)
    type(accrued_interest_result) :: out
    if (present(anniversary_dates)) then
      out = accrued_interest(start_date, end_date, coupon_rate_percent, dcc, &
        redemption_value, coupon_frequency, maturity, next_coupon_year, &
        end_of_month, anniversary_dates)
    else
      out = accrued_interest(start_date, end_date, coupon_rate_percent, dcc, &
        redemption_value, coupon_frequency, maturity, next_coupon_year, end_of_month)
    end if
  end function AccrInt

  function AnnivDates(terms) result(schedule)
    type(bond_terms), intent(in) :: terms
    type(bond_schedule) :: schedule
    call build_bond_schedule(terms, schedule)
  end function AnnivDates

  function DP(clean_price, settlement_date, terms, schedule) result(out)
    real(real_kind), intent(in) :: clean_price
    type(date_type), intent(in) :: settlement_date
    type(bond_terms), intent(in) :: terms
    type(bond_schedule), intent(in), optional :: schedule
    type(dirty_price_result) :: out
    if (present(schedule)) then
      out = dirty_price(clean_price, settlement_date, terms, schedule)
    else
      out = dirty_price(clean_price, settlement_date, terms)
    end if
  end function DP

  function BondVal_Price(yield_percent, settlement_date, terms, &
                         simple_last_period, calculation_method, schedule) result(out)
    real(real_kind), intent(in) :: yield_percent
    type(date_type), intent(in) :: settlement_date
    type(bond_terms), intent(in) :: terms
    logical, intent(in), optional :: simple_last_period
    integer, intent(in), optional :: calculation_method
    type(bond_schedule), intent(in), optional :: schedule
    type(bond_value_result) :: out
    if (present(schedule)) then
      out = bond_price(yield_percent, settlement_date, terms, &
        simple_last_period, calculation_method, schedule)
    else
      out = bond_price(yield_percent, settlement_date, terms, &
        simple_last_period, calculation_method)
    end if
  end function BondVal_Price

  function BondVal_Yield(clean_price, settlement_date, terms, &
                         simple_last_period, precision, calculation_method, schedule) result(out)
    real(real_kind), intent(in) :: clean_price
    type(date_type), intent(in) :: settlement_date
    type(bond_terms), intent(in) :: terms
    logical, intent(in), optional :: simple_last_period
    real(real_kind), intent(in), optional :: precision
    integer, intent(in), optional :: calculation_method
    type(bond_schedule), intent(in), optional :: schedule
    type(bond_value_result) :: out
    if (present(schedule)) then
      out = bond_yield(clean_price, settlement_date, terms, &
        simple_last_period, precision, calculation_method, schedule)
    else
      out = bond_yield(clean_price, settlement_date, terms, &
        simple_last_period, precision, calculation_method)
    end if
  end function BondVal_Yield

  pure integer function leap(year) result(value)
    integer, intent(in) :: year
    value = merge(1, 0, leap_year(year))
  end function leap

  pure integer function LDM(date) result(value)
    type(date_type), intent(in) :: date
    value = merge(1, 0, is_last_day_of_month(date))
  end function LDM

  pure integer function DaysInMonth(date) result(value)
    type(date_type), intent(in) :: date
    value = days_in_month(date%year, date%month)
  end function DaysInMonth

  pure integer function DaysInYear(year) result(value)
    integer, intent(in) :: year
    value = days_in_year(year)
  end function DaysInYear

  pure integer function DayDiff(first_date, second_date) result(value)
    type(date_type), intent(in) :: first_date, second_date
    value = day_diff(first_date, second_date)
  end function DayDiff

  pure function Date_LDM(date) result(value)
    type(date_type), intent(in) :: date
    type(date_type) :: value
    value = last_day_of_month(date)
  end function Date_LDM

  pure real(real_kind) function sumC(x) result(value)
    real(real_kind), intent(in) :: x(:)
    value = sum(x)
  end function sumC

  pure integer function FirstMatch(element, vector) result(position)
    integer, intent(in) :: element
    integer, intent(in) :: vector(:)
    integer :: i
    position = size(vector) + 1
    do i = 1, size(vector)
      if (vector(i) == element) then
        position = i
        return
      end if
    end do
  end function FirstMatch

  pure integer function LeapDayInside(first_date, second_date) result(value)
    type(date_type), intent(in) :: first_date, second_date
    value = merge(1, 0, leap_day_inside(first_date, second_date))
  end function LeapDayInside

  function DIST(dcc, first_date, second_date, coupon_frequency, maturity, &
                next_coupon_year, end_of_month, anniversary_dates) result(out)
    integer, intent(in) :: dcc
    type(date_type), intent(in) :: first_date, second_date
    integer, intent(in), optional :: coupon_frequency
    type(date_type), intent(in), optional :: maturity
    integer, intent(in), optional :: next_coupon_year
    logical, intent(in), optional :: end_of_month
    type(date_type), intent(in), optional :: anniversary_dates(:)
    type(daycount_result) :: out
    if (present(anniversary_dates)) then
      out = day_count_fraction(first_date, second_date, dcc, coupon_frequency, &
        maturity, end_of_month, next_coupon_year, anniversary_dates)
    else
      out = day_count_fraction(first_date, second_date, dcc, coupon_frequency, &
        maturity, end_of_month, next_coupon_year)
    end if
  end function DIST

  function PayCalc(redemption_value, coupon_rate_decimal, dcc, start_dates, &
                   end_dates, coupon_frequency, maturity, end_of_month, &
                   anniversary_dates) result(payments)
    real(real_kind), intent(in) :: redemption_value, coupon_rate_decimal
    integer, intent(in) :: dcc
    type(date_type), intent(in) :: start_dates(:), end_dates(:)
    integer, intent(in), optional :: coupon_frequency
    type(date_type), intent(in), optional :: maturity
    logical, intent(in), optional :: end_of_month
    type(date_type), intent(in), optional :: anniversary_dates(:)
    real(real_kind), allocatable :: payments(:)
    type(daycount_result) :: dc
    integer :: i
    allocate(payments(size(start_dates)))
    if (size(end_dates) /= size(start_dates)) then
      payments = 0.0_real_kind
      return
    end if
    do i = 1, size(start_dates)
      if (present(anniversary_dates)) then
        dc = day_count_fraction(start_dates(i), end_dates(i), dcc, &
          coupon_frequency, maturity, end_of_month, end_dates(i)%year, &
          anniversary_dates)
      else
        dc = day_count_fraction(start_dates(i), end_dates(i), dcc, &
          coupon_frequency, maturity, end_of_month, end_dates(i)%year)
      end if
      payments(i) = redemption_value * coupon_rate_decimal * dc%fraction
    end do
  end function PayCalc

  pure function NumToDate(number_of_days, origin) result(value)
    integer, intent(in) :: number_of_days
    type(date_type), intent(in), optional :: origin
    type(date_type) :: value
    integer :: origin_serial
    origin_serial = 0
    if (present(origin)) origin_serial = date_to_serial(origin)
    value = serial_to_date(origin_serial + number_of_days)
  end function NumToDate

  pure function CppPrevDate(date, issue_date, reference_date, coupon_frequency, &
                            end_of_month) result(value)
    type(date_type), intent(in) :: date, issue_date, reference_date
    integer, intent(in) :: coupon_frequency
    logical, intent(in) :: end_of_month
    type(date_type) :: value, middle
    integer :: gap
    if (coupon_frequency <= 0 .or. day_diff(date, issue_date) > 0) then
      value = date_type()
      return
    end if
    gap = 365 / coupon_frequency
    middle = date_type(date%year, date%month, 15)
    value = serial_to_date(date_to_serial(middle) - gap)
    if (end_of_month) then
      value = last_day_of_month(value)
    else
      value%day = min(reference_date%day, days_in_month(value%year, value%month))
    end if
  end function CppPrevDate

  pure function CppSuccDate(date, maturity_date, reference_date, &
                            coupon_frequency, end_of_month) result(value)
    type(date_type), intent(in) :: date, maturity_date, reference_date
    integer, intent(in) :: coupon_frequency
    logical, intent(in) :: end_of_month
    type(date_type) :: value, middle
    integer :: gap
    if (coupon_frequency <= 0 .or. day_diff(maturity_date, date) > 0) then
      value = date_type()
      return
    end if
    gap = 365 / coupon_frequency
    middle = date_type(date%year, date%month, 15)
    value = serial_to_date(date_to_serial(middle) + gap)
    if (end_of_month) then
      value = last_day_of_month(value)
    else
      value%day = min(reference_date%day, days_in_month(value%year, value%month))
    end if
  end function CppSuccDate

  function NewtonRaphson(function_value, derivative_value, start_value, precision, &
                         max_iterations) result(out)
    procedure(scalar_function) :: function_value, derivative_value
    real(real_kind), intent(in) :: start_value
    real(real_kind), intent(in), optional :: precision
    integer, intent(in), optional :: max_iterations
    type(newton_result) :: out
    real(real_kind) :: tolerance, current, next_value, f, derivative
    integer :: maximum, i
    tolerance = epsilon(1.0_real_kind) ** 0.25_real_kind
    if (present(precision)) tolerance = precision
    maximum = 1000
    if (present(max_iterations)) maximum = max_iterations
    current = start_value
    do i = 0, maximum
      f = function_value(current)
      if (abs(f) <= tolerance) exit
      derivative = derivative_value(current)
      if (abs(derivative) <= tiny(1.0_real_kind)) then
        out%status = 1
        exit
      end if
      next_value = current - f / derivative
      if (.not. ieee_is_finite(next_value)) then
        out%status = 2
        exit
      end if
      current = next_value
    end do
    out%root = current
    out%iterations = i
    out%function_value = function_value(current)
    if (i > maximum) out%status = 3
  end function NewtonRaphson

  real(real_kind) function dm_MyPriceEqn(a, m, cn_tau, coupon_flow, final_flow, &
                                         w, eta, z, coupon_frequency) result(value)
    real(real_kind), intent(in) :: a, cn_tau, coupon_flow, final_flow, w, eta, z
    integer, intent(in) :: m, coupon_frequency
    integer :: j
    real(real_kind) :: factor1, summand1, inner_factor, inner_sum, last_summand
    factor1 = (-1.0_real_kind) ** m / &
      (a ** (w + real(m, real_kind)) * real(coupon_frequency, real_kind) ** m)
    summand1 = cn_tau * gamma(w + real(m, real_kind)) / gamma(w)
    inner_factor = coupon_flow * real(factorial_integer(m), real_kind) / (a - 1.0_real_kind)
    inner_sum = 0.0_real_kind
    do j = 0, m
      inner_sum = inner_sum + 1.0_real_kind / real(factorial_integer(m - j), real_kind) * &
        (a / (a - 1.0_real_kind)) ** j * &
        (gamma(w + real(m - j, real_kind)) / gamma(w) - &
         gamma(w + eta + real(m - j, real_kind)) / &
         (a ** eta * gamma(w + eta)))
    end do
    last_summand = final_flow * gamma(w + eta + z + real(m, real_kind)) / &
      (a ** (eta + z) * gamma(w + eta + z))
    value = factor1 * (summand1 + inner_factor * inner_sum + last_summand)
  end function dm_MyPriceEqn

  real(real_kind) function ModDUR(a, cn_tau, coupon_flow, final_flow, w, eta, z, &
                                  coupon_frequency, dirty_price_value) result(value)
    real(real_kind), intent(in) :: a, cn_tau, coupon_flow, final_flow, w, eta, z
    integer, intent(in) :: coupon_frequency
    real(real_kind), intent(in) :: dirty_price_value
    real(real_kind) :: factor1, summand1, inner_factor, inner_sum, last_summand
    factor1 = 1.0_real_kind / (dirty_price_value * real(coupon_frequency, real_kind) * a ** (w + 1.0_real_kind))
    summand1 = cn_tau * w
    inner_factor = coupon_flow * (a ** eta - 1.0_real_kind) / &
      (a ** eta * (a - 1.0_real_kind))
    inner_sum = w - eta / (a ** eta - 1.0_real_kind) + a / (a - 1.0_real_kind)
    last_summand = final_flow * (w + eta + z) / a ** (eta + z)
    value = factor1 * (summand1 + inner_factor * inner_sum + last_summand)
  end function ModDUR

  real(real_kind) function CONV(a, cn_tau, coupon_flow, final_flow, w, eta, z, &
                                coupon_frequency, dirty_price_value) result(value)
    real(real_kind), intent(in) :: a, cn_tau, coupon_flow, final_flow, w, eta, z
    integer, intent(in) :: coupon_frequency
    real(real_kind), intent(in) :: dirty_price_value
    real(real_kind) :: factor1, summand1, inner_factor, inner_sum, last_summand
    factor1 = 0.5_real_kind / (a ** w * dirty_price_value) * &
      (1.0_real_kind / (real(coupon_frequency, real_kind) * a)) ** 2
    summand1 = cn_tau * w * (w + 1.0_real_kind)
    inner_factor = coupon_flow * (a ** eta - 1.0_real_kind) / &
      (a ** eta * (a - 1.0_real_kind))
    inner_sum = w * (w + 1.0_real_kind) - eta * (2.0_real_kind * w + eta + 1.0_real_kind) / &
      (a ** eta - 1.0_real_kind) + 2.0_real_kind * a / (a - 1.0_real_kind) * &
      (w - eta / (a ** eta - 1.0_real_kind) + a / (a - 1.0_real_kind))
    last_summand = final_flow * (w + eta + z) * (w + eta + z + 1.0_real_kind) / a ** (eta + z)
    value = factor1 * (summand1 + inner_factor * inner_sum + last_summand)
  end function CONV

  pure integer function factorial_integer(n) result(value)
    integer, intent(in) :: n
    integer :: i
    value = 1
    do i = 2, n
      value = value * i
    end do
  end function factorial_integer

end module bondvaluation_compat
