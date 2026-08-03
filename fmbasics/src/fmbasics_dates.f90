! SPDX-License-Identifier: GPL-2.0-only
module fmbasics_dates
   use fmbasics_kinds, only : dp, FM_OK, FM_INVALID_ARGUMENT
   implicit none
   private

   integer, parameter, public :: DATE_NA = -huge(0)

   type, public :: period_t
      integer :: months = 0
      integer :: days = 0
   end type period_t

   type, public :: calendar_t
      character(len=8), allocatable :: code(:)
   contains
      procedure :: size => calendar_size
   end type calendar_t

   public :: make_date, date_from_yyyymmdd, date_to_yyyymmdd
   public :: date_year, date_month, date_day, day_of_week, days_in_month
   public :: days_period, months_period, years_period
   public :: calendar, joint_calendar, calendar_contains
   public :: is_good_day, shift_date, add_months, roll_date, year_frac
   public :: date_string

   interface calendar
      module procedure calendar_scalar
      module procedure calendar_vector
   end interface calendar

   interface year_frac
      module procedure year_frac_scalar
      module procedure year_frac_vector
   end interface year_frac

contains

   pure integer function make_date(year, month, day) result(serial)
      integer, intent(in) :: year, month, day
      integer :: y, era, yoe, doy, doe, mp
      if (month < 1 .or. month > 12 .or. day < 1 .or. day > days_in_month(year, month)) then
         serial = DATE_NA
         return
      end if
      y = year
      if (month <= 2) y = y - 1
      if (y >= 0) then
         era = y / 400
      else
         era = (y - 399) / 400
      end if
      yoe = y - era * 400
      if (month > 2) then
         mp = month - 3
      else
         mp = month + 9
      end if
      doy = (153 * mp + 2) / 5 + day - 1
      doe = yoe * 365 + yoe / 4 - yoe / 100 + doy
      serial = era * 146097 + doe - 719468
   end function make_date

   pure subroutine civil_from_days(serial, year, month, day)
      integer, intent(in) :: serial
      integer, intent(out) :: year, month, day
      integer :: z, era, doe, yoe, doy, mp
      z = serial + 719468
      if (z >= 0) then
         era = z / 146097
      else
         era = (z - 146096) / 146097
      end if
      doe = z - era * 146097
      yoe = (doe - doe / 1460 + doe / 36524 - doe / 146096) / 365
      year = yoe + era * 400
      doy = doe - (365 * yoe + yoe / 4 - yoe / 100)
      mp = (5 * doy + 2) / 153
      day = doy - (153 * mp + 2) / 5 + 1
      if (mp < 10) then
         month = mp + 3
      else
         month = mp - 9
      end if
      if (month <= 2) year = year + 1
   end subroutine civil_from_days

   pure integer function date_from_yyyymmdd(value) result(serial)
      integer, intent(in) :: value
      integer :: y, m, d
      y = value / 10000
      m = mod(value / 100, 100)
      d = mod(value, 100)
      serial = make_date(y, m, d)
   end function date_from_yyyymmdd

   pure integer function date_to_yyyymmdd(serial) result(value)
      integer, intent(in) :: serial
      integer :: y, m, d
      if (serial == DATE_NA) then
         value = DATE_NA
         return
      end if
      call civil_from_days(serial, y, m, d)
      value = 10000 * y + 100 * m + d
   end function date_to_yyyymmdd

   pure integer function date_year(serial) result(value)
      integer, intent(in) :: serial
      integer :: m, d
      if (serial == DATE_NA) then
         value = DATE_NA
      else
         call civil_from_days(serial, value, m, d)
      end if
   end function date_year

   pure integer function date_month(serial) result(value)
      integer, intent(in) :: serial
      integer :: y, d
      if (serial == DATE_NA) then
         value = DATE_NA
      else
         call civil_from_days(serial, y, value, d)
      end if
   end function date_month

   pure integer function date_day(serial) result(value)
      integer, intent(in) :: serial
      integer :: y, m
      if (serial == DATE_NA) then
         value = DATE_NA
      else
         call civil_from_days(serial, y, m, value)
      end if
   end function date_day

   pure integer function days_in_month(year, month) result(value)
      integer, intent(in) :: year, month
      logical :: leap
      integer, parameter :: mdays(12) = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
      if (month < 1 .or. month > 12) then
         value = 0
         return
      end if
      leap = mod(year, 4) == 0 .and. (mod(year, 100) /= 0 .or. mod(year, 400) == 0)
      value = mdays(month)
      if (month == 2 .and. leap) value = 29
   end function days_in_month

   pure integer function day_of_week(serial) result(value)
      integer, intent(in) :: serial
      if (serial == DATE_NA) then
         value = 0
      else
         value = modulo(serial + 3, 7) + 1
      end if
   end function day_of_week

   pure function days_period(n) result(value)
      integer, intent(in) :: n
      type(period_t) :: value
      value%days = n
   end function days_period

   pure function months_period(n) result(value)
      integer, intent(in) :: n
      type(period_t) :: value
      value%months = n
   end function months_period

   pure function years_period(n) result(value)
      integer, intent(in) :: n
      type(period_t) :: value
      value%months = 12 * n
   end function years_period

   function calendar_scalar(code) result(value)
      character(len=*), intent(in) :: code
      type(calendar_t) :: value
      allocate(value%code(1))
      value%code(1) = upper(adjustl(code))
   end function calendar_scalar

   function calendar_vector(code) result(value)
      character(len=*), intent(in) :: code(:)
      type(calendar_t) :: value
      integer :: i
      allocate(value%code(size(code)))
      do i = 1, size(code)
         value%code(i) = upper(adjustl(code(i)))
      end do
   end function calendar_vector

   function joint_calendar(a, b) result(value)
      type(calendar_t), intent(in) :: a, b
      type(calendar_t) :: value
      integer :: i, n
      character(len=8), allocatable :: tmp(:)
      allocate(tmp(a%size() + b%size()))
      n = 0
      do i = 1, a%size()
         if (.not. any(tmp(:n) == a%code(i))) then
            n = n + 1
            tmp(n) = a%code(i)
         end if
      end do
      do i = 1, b%size()
         if (.not. any(tmp(:n) == b%code(i))) then
            n = n + 1
            tmp(n) = b%code(i)
         end if
      end do
      allocate(value%code(n))
      if (n > 0) value%code = tmp(:n)
   end function joint_calendar

   pure integer function calendar_size(self) result(value)
      class(calendar_t), intent(in) :: self
      if (allocated(self%code)) then
         value = size(self%code)
      else
         value = 0
      end if
   end function calendar_size

   pure logical function calendar_contains(cal, code) result(value)
      type(calendar_t), intent(in) :: cal
      character(len=*), intent(in) :: code
      if (cal%size() == 0) then
         value = .false.
      else
         value = any(cal%code == upper(adjustl(code)))
      end if
   end function calendar_contains

   pure logical function is_good_day(serial, cal) result(value)
      integer, intent(in) :: serial
      type(calendar_t), intent(in) :: cal
      integer :: i, wd
      if (serial == DATE_NA) then
         value = .false.
         return
      end if
      wd = day_of_week(serial)
      if (wd >= 6) then
         value = .false.
         return
      end if
      value = .true.
      do i = 1, cal%size()
         if (is_holiday(serial, trim(cal%code(i)))) then
            value = .false.
            return
         end if
      end do
   end function is_good_day

   pure integer function shift_date(serial, tenor, convention, cal, end_to_end) result(value)
      integer, intent(in) :: serial
      type(period_t), intent(in) :: tenor
      character(len=*), intent(in) :: convention
      type(calendar_t), intent(in) :: cal
      logical, intent(in) :: end_to_end
      integer :: i, direction
      value = serial
      if (serial == DATE_NA) return
      if (tenor%months /= 0) value = add_months(value, tenor%months, end_to_end)
      if (tenor%days /= 0) then
         direction = merge(1, -1, tenor%days > 0)
         do i = 1, abs(tenor%days)
            do
               value = value + direction
               if (is_good_day(value, cal)) exit
            end do
         end do
      end if
      value = roll_date(value, convention, cal)
   end function shift_date

   pure integer function add_months(serial, nmonths, preserve_eom) result(value)
      integer, intent(in) :: serial, nmonths
      logical, intent(in) :: preserve_eom
      integer :: y, m, d, total, ny, nm, nd
      if (serial == DATE_NA) then
         value = DATE_NA
         return
      end if
      call civil_from_days(serial, y, m, d)
      total = y * 12 + (m - 1) + nmonths
      if (total >= 0) then
         ny = total / 12
      else
         ny = (total - 11) / 12
      end if
      nm = total - ny * 12 + 1
      if (preserve_eom .and. d == days_in_month(y, m)) then
         nd = days_in_month(ny, nm)
      else
         nd = min(d, days_in_month(ny, nm))
      end if
      value = make_date(ny, nm, nd)
   end function add_months

   pure integer function roll_date(serial, convention, cal) result(value)
      integer, intent(in) :: serial
      character(len=*), intent(in) :: convention
      type(calendar_t), intent(in) :: cal
      character(len=4) :: c
      integer :: fwd, back, month0
      value = serial
      if (serial == DATE_NA .or. is_good_day(serial, cal)) return
      c = lower(adjustl(convention))
      select case (trim(c))
      case ('p', 'preceding')
         do while (.not. is_good_day(value, cal))
            value = value - 1
         end do
      case ('mf', 'mp', 'ms')
         month0 = date_month(serial)
         fwd = serial
         do while (.not. is_good_day(fwd, cal))
            fwd = fwd + 1
         end do
         if (date_month(fwd) == month0) then
            value = fwd
         else
            back = serial
            do while (.not. is_good_day(back, cal))
               back = back - 1
            end do
            value = back
         end if
      case default
         do while (.not. is_good_day(value, cal))
            value = value + 1
         end do
      end select
   end function roll_date

   real(dp) function year_frac_scalar(d1, d2, basis, status) result(value)
      integer, intent(in) :: d1, d2
      character(len=*), intent(in) :: basis
      integer, intent(out), optional :: status
      integer :: y1, m1, dd1, y2, m2, dd2
      character(len=16) :: b
      if (d1 == DATE_NA .or. d2 == DATE_NA) then
         value = 0.0_dp
         if (present(status)) status = FM_INVALID_ARGUMENT
         return
      end if
      b = lower(adjustl(basis))
      select case (trim(b))
      case ('act/365', 'actual/365', 'act365')
         value = real(d2 - d1, dp) / 365.0_dp
      case ('act/360', 'actual/360', 'act360')
         value = real(d2 - d1, dp) / 360.0_dp
      case ('act/act', 'actual/actual')
         value = actual_actual(d1, d2)
      case ('30e/360', '30/360e', 'eurobond')
         call civil_from_days(d1, y1, m1, dd1)
         call civil_from_days(d2, y2, m2, dd2)
         dd1 = min(dd1, 30)
         dd2 = min(dd2, 30)
         value = real(360 * (y2-y1) + 30 * (m2-m1) + dd2-dd1, dp) / 360.0_dp
      case ('30/360us', '30u/360', 'bond')
         call civil_from_days(d1, y1, m1, dd1)
         call civil_from_days(d2, y2, m2, dd2)
         if (dd1 == 31 .or. (m1 == 2 .and. dd1 == days_in_month(y1, m1))) dd1 = 30
         if (dd2 == 31 .and. dd1 == 30) dd2 = 30
         value = real(360 * (y2-y1) + 30 * (m2-m1) + dd2-dd1, dp) / 360.0_dp
      case default
         value = 0.0_dp
         if (present(status)) status = FM_INVALID_ARGUMENT
         return
      end select
      if (present(status)) status = FM_OK
   end function year_frac_scalar

   function year_frac_vector(d1, d2, basis, status) result(value)
      integer, intent(in) :: d1(:), d2(:)
      character(len=*), intent(in) :: basis
      integer, intent(out), optional :: status
      real(dp), allocatable :: value(:)
      integer :: i, n, stat_i
      n = max(size(d1), size(d2))
      allocate(value(n))
      stat_i = FM_OK
      do i = 1, n
         value(i) = year_frac_scalar(d1(mod(i-1,size(d1))+1), d2(mod(i-1,size(d2))+1), basis, stat_i)
         if (stat_i /= FM_OK) exit
      end do
      if (present(status)) status = stat_i
   end function year_frac_vector

   pure real(dp) function actual_actual(d1, d2) result(value)
      integer, intent(in) :: d1, d2
      integer :: start_date, end_date, y, year_end, direction
      real(dp) :: total
      if (d1 == d2) then
         value = 0.0_dp
         return
      end if
      direction = merge(1, -1, d2 > d1)
      start_date = min(d1, d2)
      end_date = max(d1, d2)
      total = 0.0_dp
      do while (start_date < end_date)
         y = date_year(start_date)
         year_end = min(end_date, make_date(y + 1, 1, 1))
         total = total + real(year_end - start_date, dp) / &
            real(merge(366, 365, days_in_month(y, 2) == 29), dp)
         start_date = year_end
      end do
      value = real(direction, dp) * total
   end function actual_actual

   pure logical function is_holiday(serial, code) result(value)
      integer, intent(in) :: serial
      character(len=*), intent(in) :: code
      integer :: y, m, d, easter, wd
      character(len=8) :: c
      call civil_from_days(serial, y, m, d)
      c = upper(adjustl(code))
      easter = easter_sunday(y)
      value = .false.

      if (observed_fixed(serial, y, 1, 1, c == 'USNY')) value = .true.
      if (serial == easter - 2 .or. serial == easter + 1) value = .true.
      if (observed_fixed(serial, y, 12, 25, c == 'USNY')) value = .true.

      select case (trim(c))
      case ('AUSY')
         value = value .or. observed_fixed(serial, y, 1, 26, .false.)
         value = value .or. observed_fixed(serial, y, 4, 25, .false.)
         value = value .or. serial == nth_weekday(y, 6, 1, 2)
         value = value .or. serial == nth_weekday(y, 10, 1, 1)
         value = value .or. observed_fixed(serial, y, 12, 26, .false.)
      case ('GBLO')
         value = value .or. serial == nth_weekday(y, 5, 1, 1)
         value = value .or. serial == last_weekday(y, 5, 1)
         value = value .or. serial == last_weekday(y, 8, 1)
         value = value .or. observed_fixed(serial, y, 12, 26, .false.)
      case ('EUTA')
         value = value .or. observed_fixed(serial, y, 5, 1, .false.)
         value = value .or. observed_fixed(serial, y, 12, 26, .false.)
      case ('USNY')
         value = value .or. serial == nth_weekday(y, 1, 1, 3)
         value = value .or. serial == nth_weekday(y, 2, 1, 3)
         value = value .or. serial == last_weekday(y, 5, 1)
         if (y >= 2022) value = value .or. observed_fixed(serial, y, 6, 19, .true.)
         value = value .or. observed_fixed(serial, y, 7, 4, .true.)
         value = value .or. serial == nth_weekday(y, 9, 1, 1)
         value = value .or. serial == nth_weekday(y, 10, 1, 2)
         value = value .or. observed_fixed(serial, y, 11, 11, .true.)
         value = value .or. serial == nth_weekday(y, 11, 4, 4)
      case ('NZAU', 'NZWE')
         value = value .or. observed_fixed(serial, y, 2, 6, .false.)
         value = value .or. observed_fixed(serial, y, 4, 25, .false.)
         value = value .or. serial == nth_weekday(y, 6, 1, 1)
         value = value .or. serial == nth_weekday(y, 10, 1, 4)
         value = value .or. observed_fixed(serial, y, 12, 26, .false.)
      case ('CHZH', 'HKHK', 'NOOS', 'JPTO')
         ! Common fixed holidays only; local holiday sets are intentionally partial.
      case default
         ! Unknown calendars are weekend-only.
      end select

      wd = day_of_week(serial)
      if (wd >= 6) value = .false.
   end function is_holiday

   pure logical function observed_fixed(serial, year, month, day, us_style) result(value)
      integer, intent(in) :: serial, year, month, day
      logical, intent(in) :: us_style
      integer :: holiday, wd, observed
      holiday = make_date(year, month, day)
      wd = day_of_week(holiday)
      observed = holiday
      if (wd == 6) then
         observed = merge(holiday - 1, holiday + 2, us_style)
      else if (wd == 7) then
         observed = holiday + 1
      end if
      value = serial == holiday .or. serial == observed
   end function observed_fixed

   pure integer function nth_weekday(year, month, weekday, nth) result(serial)
      integer, intent(in) :: year, month, weekday, nth
      integer :: first
      first = make_date(year, month, 1)
      serial = first + modulo(weekday - day_of_week(first), 7) + 7 * (nth - 1)
   end function nth_weekday

   pure integer function last_weekday(year, month, weekday) result(serial)
      integer, intent(in) :: year, month, weekday
      integer :: last
      last = make_date(year, month, days_in_month(year, month))
      serial = last - modulo(day_of_week(last) - weekday, 7)
   end function last_weekday

   pure integer function easter_sunday(year) result(serial)
      integer, intent(in) :: year
      integer :: a, b, c, d, e, f, g, h, i, k, l, m, month, day
      a = mod(year, 19)
      b = year / 100
      c = mod(year, 100)
      d = b / 4
      e = mod(b, 4)
      f = (b + 8) / 25
      g = (b - f + 1) / 3
      h = mod(19 * a + b - d - g + 15, 30)
      i = c / 4
      k = mod(c, 4)
      l = mod(32 + 2 * e + 2 * i - h - k, 7)
      m = (a + 11 * h + 22 * l) / 451
      month = (h + l - 7 * m + 114) / 31
      day = mod(h + l - 7 * m + 114, 31) + 1
      serial = make_date(year, month, day)
   end function easter_sunday

   pure function upper(text) result(out)
      character(len=*), intent(in) :: text
      character(len=len(text)) :: out
      integer :: i, k
      do i = 1, len(text)
         k = iachar(text(i:i))
         if (k >= iachar('a') .and. k <= iachar('z')) then
            out(i:i) = achar(k - 32)
         else
            out(i:i) = text(i:i)
         end if
      end do
   end function upper

   pure function lower(text) result(out)
      character(len=*), intent(in) :: text
      character(len=len(text)) :: out
      integer :: i, k
      do i = 1, len(text)
         k = iachar(text(i:i))
         if (k >= iachar('A') .and. k <= iachar('Z')) then
            out(i:i) = achar(k + 32)
         else
            out(i:i) = text(i:i)
         end if
      end do
   end function lower

   function date_string(serial) result(text)
      integer, intent(in) :: serial
      character(len=10) :: text
      integer :: y, m, d
      if (serial == DATE_NA) then
         text = 'NA        '
      else
         call civil_from_days(serial, y, m, d)
         write(text, '(i4.4,"-",i2.2,"-",i2.2)') y, m, d
      end if
   end function date_string

end module fmbasics_dates
