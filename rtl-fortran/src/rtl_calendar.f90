! SPDX-License-Identifier: MIT
! Copyright (c) 2020 RTL Authors
module rtl_calendar
  use rtl_kinds, only: dp
  use rtl_types, only: commodity_weight_result
  implicit none
  private

  public :: ymd_to_serial, serial_to_ymd, day_of_week, is_business_day
  public :: add_months, month_end, business_days_between, swap_fut_weight, swapFutWeight

interface swapFutWeight
    module procedure swap_fut_weight
  end interface swapFutWeight

contains

  pure integer function ymd_to_serial(year, month, day) result(serial)
    integer, intent(in) :: year, month, day
    integer :: y, m
    y = year
    m = month
    if (m <= 2) then
      y = y - 1
      m = m + 12
    end if
    serial = 365 * y + y / 4 - y / 100 + y / 400 + &
      (153 * (m - 3) + 2) / 5 + day - 1
  end function ymd_to_serial

  pure subroutine serial_to_ymd(serial, year, month, day)
    integer, intent(in) :: serial
    integer, intent(out) :: year, month, day
    integer :: z, era, doe, yoe, doy, mp
    z = serial
    era = floor_div(z, 146097)
    doe = z - era * 146097
    yoe = (doe - doe / 1460 + doe / 36524 - doe / 146096) / 365
    year = yoe + era * 400
    doy = doe - (365 * yoe + yoe / 4 - yoe / 100)
    mp = (5 * doy + 2) / 153
    day = doy - (153 * mp + 2) / 5 + 1
    month = mp + 3
    if (month > 12) then
      month = month - 12
      year = year + 1
    end if
  contains
    pure integer function floor_div(a, b) result(q)
      integer, intent(in) :: a, b
      q = a / b
      if (a < 0 .and. modulo(a, b) /= 0) q = q - 1
    end function floor_div
  end subroutine serial_to_ymd

  pure integer function day_of_week(serial) result(value)
    integer, intent(in) :: serial
    ! 1=Monday, ..., 7=Sunday. 1970-01-01 was Thursday.
    integer, parameter :: epoch = 719468
    value = modulo(serial - epoch + 3, 7) + 1
  end function day_of_week

  pure logical function is_business_day(serial, holidays) result(value)
    integer, intent(in) :: serial
    integer, intent(in), optional :: holidays(:)
    integer :: weekday
    weekday = day_of_week(serial)
    value = weekday <= 5
    if (value .and. present(holidays)) value = .not. any(holidays == serial)
  end function is_business_day

  pure integer function days_in_month(year, month) result(value)
    integer, intent(in) :: year, month
    integer, parameter :: month_days(12) = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    value = month_days(month)
    if (month == 2) then
      if (mod(year, 400) == 0 .or. (mod(year, 4) == 0 .and. mod(year, 100) /= 0)) value = 29
    end if
  end function days_in_month

  pure integer function add_months(serial, number_months) result(value)
    integer, intent(in) :: serial, number_months
    integer :: year, month, day, total_months, new_year, new_month, new_day
    call serial_to_ymd(serial, year, month, day)
    total_months = year * 12 + month - 1 + number_months
    new_year = total_months / 12
    new_month = modulo(total_months, 12) + 1
    new_day = min(day, days_in_month(new_year, new_month))
    value = ymd_to_serial(new_year, new_month, new_day)
  end function add_months

  pure integer function month_end(serial) result(value)
    integer, intent(in) :: serial
    integer :: year, month, day
    call serial_to_ymd(serial, year, month, day)
    value = ymd_to_serial(year, month, days_in_month(year, month))
  end function month_end

  pure integer function business_days_between(start_date, end_date, holidays) result(count)
    integer, intent(in) :: start_date, end_date
    integer, intent(in), optional :: holidays(:)
    integer :: date_value
    count = 0
    do date_value = start_date, end_date
      if (present(holidays)) then
        if (is_business_day(date_value, holidays)) count = count + 1
      else
        if (is_business_day(date_value)) count = count + 1
      end if
    end do
  end function business_days_between

  function swap_fut_weight(month_start, expiry_date, holidays) result(output)
    integer, intent(in) :: month_start, expiry_date
    integer, intent(in), optional :: holidays(:)
    type(commodity_weight_result) :: output
    integer :: final_date, current_date, total_days

    final_date = month_end(month_start)
    total_days = 0
    do current_date = month_start, final_date
      if (present(holidays)) then
        if (.not. is_business_day(current_date, holidays)) cycle
      else
        if (.not. is_business_day(current_date)) cycle
      end if
      total_days = total_days + 1
      if (current_date <= expiry_date) then
        output%days_first = output%days_first + 1
      else
        output%days_second = output%days_second + 1
      end if
    end do
    if (total_days <= 0) then
      output%status%ok = .false.
      output%status%message = "month contains no business days"
      return
    end if
    output%first_weight = real(output%days_first, dp) / real(total_days, dp)
  end function swap_fut_weight

end module rtl_calendar
