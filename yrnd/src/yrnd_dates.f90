! SPDX-License-Identifier: GPL-3.0-only
module yrnd_dates
   use yrnd_kinds, only : dp
   implicit none
   private

   integer, parameter, public :: dc_act_act = 1
   integer, parameter, public :: dc_act_360 = 2
   integer, parameter, public :: dc_act_365 = 3
   integer, parameter, public :: dc_30_360 = 4

   type, public :: date_t
      integer :: year = 1970
      integer :: month = 1
      integer :: day = 1
   end type date_t

   interface operator(<)
      module procedure date_less
   end interface
   interface operator(<=)
      module procedure date_less_equal
   end interface
   interface operator(>)
      module procedure date_greater
   end interface
   interface operator(>=)
      module procedure date_greater_equal
   end interface
   interface operator(==)
      module procedure date_equal
   end interface

   public :: serial_day, date_from_serial, add_days, add_months, days_between
   public :: days_in_month, is_leap_year, year_fraction
   public :: operator(<), operator(<=), operator(>), operator(>=), operator(==)

contains

   pure logical function is_leap_year(year) result(leap)
      integer, intent(in) :: year
      leap = mod(year, 4) == 0 .and. (mod(year, 100) /= 0 .or. mod(year, 400) == 0)
   end function is_leap_year

   pure integer function days_in_month(year, month) result(n)
      integer, intent(in) :: year, month
      integer, parameter :: mdays(12) = [31,28,31,30,31,30,31,31,30,31,30,31]
      if (month < 1 .or. month > 12) then
         n = 0
      else
         n = mdays(month)
         if (month == 2 .and. is_leap_year(year)) n = 29
      end if
   end function days_in_month

   pure integer function serial_day(date) result(z)
      type(date_t), intent(in) :: date
      integer :: y, m, era, yoe, doy, doe
      y = date%year
      m = date%month
      if (m <= 2) y = y - 1
      if (y >= 0) then
         era = y / 400
      else
         era = (y - 399) / 400
      end if
      yoe = y - era * 400
      if (m > 2) then
         doy = (153 * (m - 3) + 2) / 5 + date%day - 1
      else
         doy = (153 * (m + 9) + 2) / 5 + date%day - 1
      end if
      doe = yoe * 365 + yoe / 4 - yoe / 100 + doy
      z = era * 146097 + doe - 719468
   end function serial_day

   pure function date_from_serial(z) result(date)
      integer, intent(in) :: z
      type(date_t) :: date
      integer :: zz, era, doe, yoe, y, doy, mp
      zz = z + 719468
      if (zz >= 0) then
         era = zz / 146097
      else
         era = (zz - 146096) / 146097
      end if
      doe = zz - era * 146097
      yoe = (doe - doe / 1460 + doe / 36524 - doe / 146096) / 365
      y = yoe + era * 400
      doy = doe - (365 * yoe + yoe / 4 - yoe / 100)
      mp = (5 * doy + 2) / 153
      date%day = doy - (153 * mp + 2) / 5 + 1
      if (mp < 10) then
         date%month = mp + 3
      else
         date%month = mp - 9
      end if
      if (date%month <= 2) y = y + 1
      date%year = y
   end function date_from_serial

   pure function add_days(date, n) result(out)
      type(date_t), intent(in) :: date
      integer, intent(in) :: n
      type(date_t) :: out
      out = date_from_serial(serial_day(date) + n)
   end function add_days

   pure function add_months(date, nmonths) result(out)
      type(date_t), intent(in) :: date
      integer, intent(in) :: nmonths
      type(date_t) :: out
      integer :: total, y, m
      total = date%year * 12 + date%month - 1 + nmonths
      if (total >= 0) then
         y = total / 12
      else
         y = (total - 11) / 12
      end if
      m = total - 12 * y + 1
      out%year = y
      out%month = m
      out%day = min(date%day, days_in_month(y, m))
   end function add_months

   pure integer function days_between(date1, date2) result(n)
      type(date_t), intent(in) :: date1, date2
      n = serial_day(date2) - serial_day(date1)
   end function days_between

   pure real(dp) function year_fraction(date1, date2, convention) result(value)
      type(date_t), intent(in) :: date1, date2
      integer, intent(in) :: convention
      integer :: d1, d2, m1, m2
      select case (convention)
      case (dc_act_act)
         value = real(days_between(date1, date2), dp) / &
            real(days_between(date_t(date1%year, 1, 1), date_t(date2%year + 1, 1, 1)), dp)
      case (dc_act_360)
         value = real(days_between(date1, date2), dp) / 360.0_dp
      case (dc_act_365)
         value = real(days_between(date1, date2), dp) / 365.0_dp
      case default
         d1 = min(date1%day, 30)
         d2 = min(date2%day, 30)
         m1 = date1%month
         m2 = date2%month
         value = real(360 * (date2%year - date1%year) + 30 * (m2 - m1) + d2 - d1, dp) / 360.0_dp
      end select
   end function year_fraction

   pure logical function date_less(a, b)
      type(date_t), intent(in) :: a, b
      date_less = serial_day(a) < serial_day(b)
   end function date_less

   pure logical function date_less_equal(a, b)
      type(date_t), intent(in) :: a, b
      date_less_equal = serial_day(a) <= serial_day(b)
   end function date_less_equal

   pure logical function date_greater(a, b)
      type(date_t), intent(in) :: a, b
      date_greater = serial_day(a) > serial_day(b)
   end function date_greater

   pure logical function date_greater_equal(a, b)
      type(date_t), intent(in) :: a, b
      date_greater_equal = serial_day(a) >= serial_day(b)
   end function date_greater_equal

   pure logical function date_equal(a, b)
      type(date_t), intent(in) :: a, b
      date_equal = serial_day(a) == serial_day(b)
   end function date_equal

end module yrnd_dates
