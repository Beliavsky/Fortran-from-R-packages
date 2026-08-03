! SPDX-License-Identifier: GPL-3.0-or-later
module qbc_dates
   use qbc_kinds, only : dp
   use qbc_status, only : qbc_success, qbc_invalid_argument
   implicit none
   private

   type, public :: qbc_date
      integer :: year = 1970
      integer :: month = 1
      integer :: day = 1
   end type qbc_date

   public :: make_date, parse_date, date_string, valid_date
   public :: days_between, add_days, add_months, add_years
   public :: is_leap_year, days_in_month, day_of_week, is_business_day
   public :: adjust_business_day, discount_time, year_fraction
   public :: operator(==), operator(/=), operator(<), operator(<=), operator(>), operator(>=)

   interface operator(==); module procedure date_eq; end interface
   interface operator(/=); module procedure date_ne; end interface
   interface operator(<); module procedure date_lt; end interface
   interface operator(<=); module procedure date_le; end interface
   interface operator(>); module procedure date_gt; end interface
   interface operator(>=); module procedure date_ge; end interface

contains

   pure function make_date(year, month, day) result(date)
      integer, intent(in) :: year, month, day
      type(qbc_date) :: date
      date%year = year
      date%month = month
      date%day = day
   end function make_date

   pure logical function is_leap_year(year)
      integer, intent(in) :: year
      is_leap_year = mod(year, 4) == 0 .and. (mod(year, 100) /= 0 .or. mod(year, 400) == 0)
   end function is_leap_year

   pure integer function days_in_month(year, month)
      integer, intent(in) :: year, month
      integer, parameter :: mdays(12) = [31,28,31,30,31,30,31,31,30,31,30,31]
      if (month < 1 .or. month > 12) then
         days_in_month = 0
      else
         days_in_month = mdays(month)
         if (month == 2 .and. is_leap_year(year)) days_in_month = 29
      end if
   end function days_in_month

   pure logical function valid_date(date)
      type(qbc_date), intent(in) :: date
      valid_date = date%month >= 1 .and. date%month <= 12
      if (valid_date) valid_date = date%day >= 1 .and. date%day <= days_in_month(date%year, date%month)
   end function valid_date

   subroutine parse_date(text, date, status)
      character(len=*), intent(in) :: text
      type(qbc_date), intent(out) :: date
      integer, intent(out), optional :: status
      integer :: ios, st
      character(len=:), allocatable :: s
      st = qbc_success
      s = adjustl(trim(text))
      if (len_trim(s) /= 10 .or. s(5:5) /= '-' .or. s(8:8) /= '-') then
         st = qbc_invalid_argument
         date = qbc_date()
      else
         read(s(1:4), *, iostat=ios) date%year
         if (ios == 0) read(s(6:7), *, iostat=ios) date%month
         if (ios == 0) read(s(9:10), *, iostat=ios) date%day
         if (ios /= 0 .or. .not. valid_date(date)) then
            st = qbc_invalid_argument
            date = qbc_date()
         end if
      end if
      if (present(status)) status = st
   end subroutine parse_date

   pure function date_string(date) result(text)
      type(qbc_date), intent(in) :: date
      character(len=10) :: text
      write(text, '(i4.4,"-",i2.2,"-",i2.2)') date%year, date%month, date%day
   end function date_string

   pure integer function serial_day(date) result(z)
      type(qbc_date), intent(in) :: date
      integer :: y, m, era, yoe, doy, doe
      y = date%year
      m = date%month
      if (m <= 2) y = y - 1
      era = floor(real(y, dp) / 400.0_dp)
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
      type(qbc_date) :: date
      integer :: zz, era, doe, yoe, y, doy, mp
      zz = z + 719468
      era = floor(real(zz, dp) / 146097.0_dp)
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

   pure integer function days_between(initial, final)
      type(qbc_date), intent(in) :: initial, final
      days_between = serial_day(final) - serial_day(initial)
   end function days_between

   pure function add_days(date, n) result(out)
      type(qbc_date), intent(in) :: date
      integer, intent(in) :: n
      type(qbc_date) :: out
      out = date_from_serial(serial_day(date) + n)
   end function add_days

   pure function add_months(date, nmonths) result(out)
      type(qbc_date), intent(in) :: date
      integer, intent(in) :: nmonths
      type(qbc_date) :: out
      integer :: total, y, m
      total = date%year * 12 + date%month - 1 + nmonths
      y = floor(real(total, dp) / 12.0_dp)
      m = total - y * 12 + 1
      out%year = y
      out%month = m
      out%day = min(date%day, days_in_month(y, m))
   end function add_months

   pure function add_years(date, nyears) result(out)
      type(qbc_date), intent(in) :: date
      integer, intent(in) :: nyears
      type(qbc_date) :: out
      out = add_months(date, 12 * nyears)
   end function add_years

   pure integer function day_of_week(date)
      type(qbc_date), intent(in) :: date
      ! Monday=1, ..., Sunday=7. 1970-01-01 was Thursday.
      day_of_week = modulo(serial_day(date) + 3, 7) + 1
   end function day_of_week

   pure logical function is_business_day(date)
      type(qbc_date), intent(in) :: date
      integer :: dow
      dow = day_of_week(date)
      is_business_day = dow <= 5
   end function is_business_day

   pure function next_business_day(date) result(out)
      type(qbc_date), intent(in) :: date
      type(qbc_date) :: out
      out = date
      do while (.not. is_business_day(out))
         out = add_days(out, 1)
      end do
   end function next_business_day

   pure function previous_business_day(date) result(out)
      type(qbc_date), intent(in) :: date
      type(qbc_date) :: out
      out = date
      do while (.not. is_business_day(out))
         out = add_days(out, -1)
      end do
   end function previous_business_day

   pure function adjust_business_day(date, convention) result(out)
      type(qbc_date), intent(in) :: date
      character(len=*), intent(in), optional :: convention
      type(qbc_date) :: out, trial
      character(len=2) :: conv
      conv = 'F '
      if (present(convention)) conv = adjustl(convention)
      select case (trim(conv))
      case ('F')
         out = next_business_day(date)
      case ('B')
         out = previous_business_day(date)
      case ('MF')
         trial = next_business_day(date)
         if (trial%month /= date%month) then
            out = previous_business_day(date)
         else
            out = trial
         end if
      case ('MB')
         trial = previous_business_day(date)
         if (trial%month /= date%month) then
            out = next_business_day(date)
         else
            out = trial
         end if
      case default
         out = date
      end select
   end function adjust_business_day

   real(dp) function discount_time(initial, final, status) result(value)
      type(qbc_date), intent(in) :: initial, final
      integer, intent(out), optional :: status
      integer :: years, st, denom
      type(qbc_date) :: anniversary, next_anniversary
      st = qbc_success
      value = 0.0_dp
      if (final < initial) then
         st = qbc_invalid_argument
      else
         years = max(0, final%year - initial%year)
         anniversary = add_years(initial, years)
         if (anniversary > final) then
            years = years - 1
            anniversary = add_years(initial, years)
         end if
         next_anniversary = add_years(anniversary, 1)
         denom = days_between(anniversary, next_anniversary)
         value = real(years, dp)
         if (denom > 0) value = value + real(days_between(anniversary, final), dp) / real(denom, dp)
      end if
      if (present(status)) status = st
   end function discount_time

   pure logical function interval_contains_feb29(initial, final)
      type(qbc_date), intent(in) :: initial, final
      integer :: y
      type(qbc_date) :: leapday
      interval_contains_feb29 = .false.
      do y = initial%year, final%year
         if (is_leap_year(y)) then
            leapday = make_date(y, 2, 29)
            if (leapday > initial .and. leapday <= final) then
               interval_contains_feb29 = .true.
               return
            end if
         end if
      end do
   end function interval_contains_feb29

   pure integer function count_feb29(initial, final)
      type(qbc_date), intent(in) :: initial, final
      integer :: y
      type(qbc_date) :: leapday
      count_feb29 = 0
      do y = initial%year, final%year
         if (is_leap_year(y)) then
            leapday = make_date(y, 2, 29)
            if (leapday > initial .and. leapday <= final) count_feb29 = count_feb29 + 1
         end if
      end do
   end function count_feb29

   real(dp) function year_fraction(initial, final, convention, status) result(value)
      type(qbc_date), intent(in) :: initial, final
      character(len=*), intent(in), optional :: convention
      integer, intent(out), optional :: status
      character(len=16) :: conv
      integer :: ndays, d1, d2, years, y, st
      type(qbc_date) :: a, b, ystart, yend
      real(dp) :: acc
      st = qbc_success
      value = 0.0_dp
      if (final < initial) then
         st = qbc_invalid_argument
         if (present(status)) status = st
         return
      end if
      conv = 'ACT/360'
      if (present(convention)) conv = adjustl(convention)
      ndays = days_between(initial, final)
      select case (trim(conv))
      case ('30/360', '30E/360')
         d1 = min(initial%day, 30)
         d2 = min(final%day, 30)
         value = real(360 * (final%year - initial%year) + 30 * (final%month - initial%month) + d2 - d1, dp) / 360.0_dp
      case ('ACT/365')
         value = real(ndays, dp) / 365.0_dp
      case ('ACT/360')
         value = real(ndays, dp) / 360.0_dp
      case ('ACT/365L')
         if (interval_contains_feb29(initial, final)) then
            value = real(ndays, dp) / 366.0_dp
         else
            value = real(ndays, dp) / 365.0_dp
         end if
      case ('NL/365')
         value = real(ndays - count_feb29(initial, final), dp) / 365.0_dp
      case ('ACT/ACT-AFB')
         years = 0
         a = initial
         do
            b = add_years(a, 1)
            if (b <= final) then
               years = years + 1
               a = b
            else
               exit
            end if
         end do
         if (interval_contains_feb29(a, final)) then
            value = real(years, dp) + real(days_between(a, final), dp) / 366.0_dp
         else
            value = real(years, dp) + real(days_between(a, final), dp) / 365.0_dp
         end if
      case ('ACT/ACT-ISDA')
         if (initial == final) then
            value = 0.0_dp
         else
            acc = 0.0_dp
            do y = initial%year, final%year
               if (y == initial%year) then
                  ystart = initial
               else
                  ystart = make_date(y, 1, 1)
               end if
               if (y == final%year) then
                  yend = final
               else
                  yend = make_date(y + 1, 1, 1)
               end if
               if (is_leap_year(y)) then
                  acc = acc + real(days_between(ystart, yend), dp) / 366.0_dp
               else
                  acc = acc + real(days_between(ystart, yend), dp) / 365.0_dp
               end if
            end do
            value = acc
         end if
      case default
         value = discount_time(initial, final, st)
      end select
      if (present(status)) status = st
   end function year_fraction

   pure logical function date_eq(a, b)
      type(qbc_date), intent(in) :: a, b
      date_eq = serial_day(a) == serial_day(b)
   end function date_eq
   pure logical function date_ne(a, b)
      type(qbc_date), intent(in) :: a, b
      date_ne = .not. date_eq(a, b)
   end function date_ne
   pure logical function date_lt(a, b)
      type(qbc_date), intent(in) :: a, b
      date_lt = serial_day(a) < serial_day(b)
   end function date_lt
   pure logical function date_le(a, b)
      type(qbc_date), intent(in) :: a, b
      date_le = serial_day(a) <= serial_day(b)
   end function date_le
   pure logical function date_gt(a, b)
      type(qbc_date), intent(in) :: a, b
      date_gt = serial_day(a) > serial_day(b)
   end function date_gt
   pure logical function date_ge(a, b)
      type(qbc_date), intent(in) :: a, b
      date_ge = serial_day(a) >= serial_day(b)
   end function date_ge

end module qbc_dates
