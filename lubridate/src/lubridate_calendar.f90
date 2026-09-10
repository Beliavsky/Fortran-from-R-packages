! SPDX-License-Identifier: MIT
! SPDX-FileComment: Gregorian calendar algorithms adapted from Beliavsky/DataFrame under MIT.
module lubridate_calendar
   !! Implements civil construction, serial conversion, accessors, and comparisons.
   use, intrinsic :: iso_fortran_env, only: int64
   use lubridate_types, only: date_type, datetime_type, dp
   implicit none
   private

   public :: as_date, as_datetime, date, day, days_in_month, decimal_date
   public :: date_decimal, hour, isoweek, isoyear, leap_year, make_date, make_datetime
   public :: mday, minute, month, origin, quarter, qday, second, semester
   public :: week, wday, yday, year, epiweek, epiyear
   public :: date_from_day_number, date_to_day_number, datetime_from_epoch, datetime_to_epoch
   public :: operator(+), operator(-), operator(==), operator(/=), operator(<), operator(<=)
   public :: operator(>), operator(>=)

   interface as_date
      module procedure as_date_date, as_date_datetime
   end interface as_date

   interface as_datetime
      module procedure as_datetime_date, as_datetime_datetime
   end interface as_datetime

   interface year
      module procedure year_date, year_datetime
   end interface year

   interface month
      module procedure month_date, month_datetime
   end interface month

   interface day
      module procedure day_date, day_datetime
   end interface day

   interface mday
      module procedure day_date, day_datetime
   end interface mday

   interface yday
      module procedure yday_date, yday_datetime
   end interface yday

   interface qday
      module procedure qday_date, qday_datetime
   end interface qday

   interface wday
      module procedure wday_date, wday_datetime
   end interface wday

   interface week
      module procedure week_date, week_datetime
   end interface week

   interface isoweek
      module procedure isoweek_date, isoweek_datetime
   end interface isoweek

   interface isoyear
      module procedure isoyear_date, isoyear_datetime
   end interface isoyear

   interface epiweek
      module procedure epiweek_date, epiweek_datetime
   end interface epiweek

   interface epiyear
      module procedure epiyear_date, epiyear_datetime
   end interface epiyear

   interface quarter
      module procedure quarter_date, quarter_datetime
   end interface quarter

   interface semester
      module procedure semester_date, semester_datetime
   end interface semester

   interface days_in_month
      module procedure days_in_month_date, days_in_month_datetime, days_in_month_components
   end interface days_in_month

   interface decimal_date
      module procedure decimal_date_date, decimal_date_datetime
   end interface decimal_date

   interface operator(+)
      module procedure date_add_days, days_add_date
   end interface operator(+)

   interface operator(-)
      module procedure date_subtract_days, date_difference
   end interface operator(-)

   interface operator(==)
      module procedure date_equal, datetime_equal
   end interface operator(==)

   interface operator(/=)
      module procedure date_not_equal, datetime_not_equal
   end interface operator(/=)

   interface operator(<)
      module procedure date_less, datetime_less
   end interface operator(<)

   interface operator(<=)
      module procedure date_less_equal, datetime_less_equal
   end interface operator(<=)

   interface operator(>)
      module procedure date_greater, datetime_greater
   end interface operator(>)

   interface operator(>=)
      module procedure date_greater_equal, datetime_greater_equal
   end interface operator(>=)

contains

   pure elemental function make_date(year, month, day) result(value)
      !! Constructs a validated Gregorian date, returning an invalid zero date on failure.
      integer, intent(in) :: year  !! Gregorian year.
      integer, intent(in) :: month !! Month number.
      integer, intent(in) :: day   !! Day of month.
      type(date_type) :: value

      value = date_type(year, month, day)
      if (.not. value%valid()) value = date_type()
   end function make_date

   pure elemental function make_datetime(year, month, day, hour, minute, second, offset_minutes) result(value)
      !! Constructs a validated fixed-offset date-time.
      integer, intent(in) :: year, month, day !! Civil date components.
      integer, intent(in), optional :: hour, minute !! Civil clock components.
      real(dp), intent(in), optional :: second !! Seconds including fraction.
      integer, intent(in), optional :: offset_minutes !! Minutes east of UTC.
      type(datetime_type) :: value

      value%date = make_date(year, month, day)
      if (present(hour)) value%hour = hour
      if (present(minute)) value%minute = minute
      if (present(second)) value%second = second
      if (present(offset_minutes)) value%offset_minutes = offset_minutes
      if (.not. value%valid()) value = datetime_type()
   end function make_datetime

   pure elemental function origin() result(value)
      !! Returns the Unix epoch date 1970-01-01.
      type(date_type) :: value

      value = date_type(1970, 1, 1)
   end function origin

   pure elemental integer(int64) function date_to_day_number(value) result(number)
      !! Converts a valid date to days since 1970-01-01.
      type(date_type), intent(in) :: value !! Date to encode.
      integer(int64) :: day_of_era, day_of_year_zero, era, month_prime, year, year_of_era

      if (.not. value%valid()) then
         number = -huge(0_int64)
         return
      end if
      year = int(value%year, int64)
      if (value%month <= 2) year = year - 1_int64
      era = floor_divide(year, 400_int64)
      year_of_era = year - era * 400_int64
      if (value%month > 2) then
         month_prime = int(value%month - 3, int64)
      else
         month_prime = int(value%month + 9, int64)
      end if
      day_of_year_zero = (153_int64 * month_prime + 2_int64) / 5_int64 + int(value%day - 1, int64)
      day_of_era = year_of_era * 365_int64 + year_of_era / 4_int64 - year_of_era / 100_int64 + day_of_year_zero
      number = era * 146097_int64 + day_of_era - 719468_int64
   end function date_to_day_number

   pure elemental function date_from_day_number(number) result(value)
      !! Converts days since 1970-01-01 to a proleptic-Gregorian date.
      integer(int64), intent(in) :: number !! Serial day number.
      type(date_type) :: value
      integer(int64) :: adjusted, day_of_era, day_of_year_zero, era
      integer(int64) :: month, month_prime, year, year_of_era

      adjusted = number + 719468_int64
      era = floor_divide(adjusted, 146097_int64)
      day_of_era = adjusted - era * 146097_int64
      year_of_era = (day_of_era - day_of_era / 1460_int64 + day_of_era / 36524_int64 - &
         day_of_era / 146096_int64) / 365_int64
      year = year_of_era + era * 400_int64
      day_of_year_zero = day_of_era - (365_int64 * year_of_era + year_of_era / 4_int64 - year_of_era / 100_int64)
      month_prime = (5_int64 * day_of_year_zero + 2_int64) / 153_int64
      value%day = int(day_of_year_zero - (153_int64 * month_prime + 2_int64) / 5_int64 + 1_int64)
      if (month_prime < 10_int64) then
         month = month_prime + 3_int64
      else
         month = month_prime - 9_int64
      end if
      if (month <= 2_int64) year = year + 1_int64
      value%year = int(year)
      value%month = int(month)
   end function date_from_day_number

   pure elemental real(dp) function datetime_to_epoch(value) result(seconds_since_epoch)
      !! Converts a valid fixed-offset date-time to a UTC Unix timestamp.
      type(datetime_type), intent(in) :: value !! Date-time to encode.

      if (.not. value%valid()) then
         seconds_since_epoch = -huge(0.0_dp)
      else
         seconds_since_epoch = real(date_to_day_number(value%date), dp) * 86400.0_dp
         seconds_since_epoch = seconds_since_epoch + real(value%hour * 3600 + value%minute * 60, dp) + value%second
         seconds_since_epoch = seconds_since_epoch - real(value%offset_minutes * 60, dp)
      end if
   end function datetime_to_epoch

   pure elemental function datetime_from_epoch(seconds_since_epoch, offset_minutes) result(value)
      !! Converts a UTC Unix timestamp to civil components at a fixed offset.
      real(dp), intent(in) :: seconds_since_epoch !! UTC seconds since 1970-01-01.
      integer, intent(in), optional :: offset_minutes !! Minutes east of UTC.
      type(datetime_type) :: value
      real(dp) :: local_seconds, seconds_of_day
      integer(int64) :: day_number

      if (present(offset_minutes)) value%offset_minutes = offset_minutes
      local_seconds = seconds_since_epoch + real(value%offset_minutes * 60, dp)
      day_number = int(floor(local_seconds / 86400.0_dp), int64)
      seconds_of_day = local_seconds - real(day_number, dp) * 86400.0_dp
      value%date = date_from_day_number(day_number)
      value%hour = int(seconds_of_day / 3600.0_dp)
      seconds_of_day = seconds_of_day - real(value%hour * 3600, dp)
      value%minute = int(seconds_of_day / 60.0_dp)
      value%second = seconds_of_day - real(value%minute * 60, dp)
   end function datetime_from_epoch

   pure elemental function as_date_date(value) result(output)
      !! Returns an existing date unchanged.
      type(date_type), intent(in) :: value !! Date to convert.
      type(date_type) :: output

      output = value
   end function as_date_date

   pure elemental function as_date_datetime(value) result(output)
      !! Extracts the displayed civil date from a date-time.
      type(datetime_type), intent(in) :: value !! Date-time to convert.
      type(date_type) :: output

      output = value%date
   end function as_date_datetime

   pure elemental function as_datetime_date(value, offset_minutes) result(output)
      !! Converts a date to midnight at a fixed offset.
      type(date_type), intent(in) :: value !! Date to convert.
      integer, intent(in), optional :: offset_minutes !! Minutes east of UTC.
      type(datetime_type) :: output

      output%date = value
      if (present(offset_minutes)) output%offset_minutes = offset_minutes
   end function as_datetime_date

   pure elemental function as_datetime_datetime(value, offset_minutes) result(output)
      !! Returns a date-time, optionally displaying the same instant at a new offset.
      type(datetime_type), intent(in) :: value !! Date-time to convert.
      integer, intent(in), optional :: offset_minutes !! New fixed offset.
      type(datetime_type) :: output

      output = value
      if (present(offset_minutes)) output = datetime_from_epoch(datetime_to_epoch(value), offset_minutes)
   end function as_datetime_datetime

   pure elemental function date(value) result(output)
      !! Extracts the civil date component of a date-time.
      type(datetime_type), intent(in) :: value !! Date-time to inspect.
      type(date_type) :: output

      output = value%date
   end function date

   pure elemental integer function year_date(value) result(component)
      !! Returns a date's year.
      type(date_type), intent(in) :: value !! Date to inspect.
      component = value%year
   end function year_date

   pure elemental integer function year_datetime(value) result(component)
      !! Returns a date-time's displayed year.
      type(datetime_type), intent(in) :: value !! Date-time to inspect.
      component = value%date%year
   end function year_datetime

   pure elemental integer function month_date(value) result(component)
      !! Returns a date's month number.
      type(date_type), intent(in) :: value !! Date to inspect.
      component = value%month
   end function month_date

   pure elemental integer function month_datetime(value) result(component)
      !! Returns a date-time's displayed month number.
      type(datetime_type), intent(in) :: value !! Date-time to inspect.
      component = value%date%month
   end function month_datetime

   pure elemental integer function day_date(value) result(component)
      !! Returns a date's day of month.
      type(date_type), intent(in) :: value !! Date to inspect.
      component = value%day
   end function day_date

   pure elemental integer function day_datetime(value) result(component)
      !! Returns a date-time's displayed day of month.
      type(datetime_type), intent(in) :: value !! Date-time to inspect.
      component = value%date%day
   end function day_datetime

   pure elemental integer function hour(value) result(component)
      !! Returns a date-time's displayed hour.
      type(datetime_type), intent(in) :: value !! Date-time to inspect.
      component = value%hour
   end function hour

   pure elemental integer function minute(value) result(component)
      !! Returns a date-time's displayed minute.
      type(datetime_type), intent(in) :: value !! Date-time to inspect.
      component = value%minute
   end function minute

   pure elemental real(dp) function second(value) result(component)
      !! Returns a date-time's displayed seconds including fraction.
      type(datetime_type), intent(in) :: value !! Date-time to inspect.
      component = value%second
   end function second

   pure elemental integer function yday_date(value) result(component)
      !! Returns the one-based day of year.
      type(date_type), intent(in) :: value !! Date to inspect.

      component = int(date_to_day_number(value) - date_to_day_number(date_type(value%year, 1, 1)) + 1_int64)
   end function yday_date

   pure elemental integer function yday_datetime(value) result(component)
      !! Returns a date-time's one-based day of year.
      type(datetime_type), intent(in) :: value !! Date-time to inspect.
      component = yday_date(value%date)
   end function yday_datetime

   pure elemental integer function qday_date(value) result(component)
      !! Returns the one-based day within the calendar quarter.
      type(date_type), intent(in) :: value !! Date to inspect.
      integer :: quarter_start_month

      quarter_start_month = 3 * ((value%month - 1) / 3) + 1
      component = int(date_to_day_number(value) - &
         date_to_day_number(date_type(value%year, quarter_start_month, 1)) + 1_int64)
   end function qday_date

   pure elemental integer function qday_datetime(value) result(component)
      !! Returns a date-time's one-based day within its quarter.
      type(datetime_type), intent(in) :: value !! Date-time to inspect.
      component = qday_date(value%date)
   end function qday_datetime

   pure elemental integer function wday_date(value, week_start) result(component)
      !! Returns weekday numbering from a caller-selected week start.
      type(date_type), intent(in) :: value !! Date to inspect.
      integer, intent(in), optional :: week_start !! R convention 1=Sunday through 7=Saturday.
      integer :: start, sunday_number

      start = 1
      if (present(week_start)) start = week_start
      sunday_number = modulo(int(date_to_day_number(value)) + 4, 7) + 1
      component = modulo(sunday_number - start, 7) + 1
   end function wday_date

   pure elemental integer function wday_datetime(value, week_start) result(component)
      !! Returns a date-time's weekday number.
      type(datetime_type), intent(in) :: value !! Date-time to inspect.
      integer, intent(in), optional :: week_start !! Week-start convention.
      component = wday_date(value%date, week_start)
   end function wday_datetime

   pure elemental integer function week_date(value) result(component)
      !! Returns lubridate's simple seven-day-block week of year.
      type(date_type), intent(in) :: value !! Date to inspect.
      component = (yday_date(value) - 1) / 7 + 1
   end function week_date

   pure elemental integer function week_datetime(value) result(component)
      !! Returns a date-time's simple week of year.
      type(datetime_type), intent(in) :: value !! Date-time to inspect.
      component = week_date(value%date)
   end function week_datetime

   pure elemental integer function isoweek_date(value) result(component)
      !! Returns the ISO-8601 week number.
      type(date_type), intent(in) :: value !! Date to inspect.
      integer(int64) :: first_monday, serial
      integer :: iso_year_value

      serial = date_to_day_number(value)
      iso_year_value = isoyear_date(value)
      first_monday = date_to_day_number(date_type(iso_year_value, 1, 4))
      first_monday = first_monday - int(wday_date(date_type(iso_year_value, 1, 4), week_start=2) - 1, int64)
      component = int((serial - first_monday) / 7_int64 + 1_int64)
   end function isoweek_date

   pure elemental integer function isoweek_datetime(value) result(component)
      !! Returns a date-time's ISO-8601 week number.
      type(datetime_type), intent(in) :: value !! Date-time to inspect.
      component = isoweek_date(value%date)
   end function isoweek_datetime

   pure elemental integer function isoyear_date(value) result(component)
      !! Returns the ISO-8601 week-numbering year.
      type(date_type), intent(in) :: value !! Date to inspect.
      type(date_type) :: thursday

      thursday = date_from_day_number(date_to_day_number(value) + &
         int(4 - wday_date(value, week_start=2), int64))
      component = thursday%year
   end function isoyear_date

   pure elemental integer function isoyear_datetime(value) result(component)
      !! Returns a date-time's ISO week-numbering year.
      type(datetime_type), intent(in) :: value !! Date-time to inspect.
      component = isoyear_date(value%date)
   end function isoyear_datetime

   pure elemental integer function epiweek_date(value) result(component)
      !! Returns the CDC epidemiological week number using Sunday-start weeks.
      type(date_type), intent(in) :: value !! Date to inspect.
      integer(int64) :: first_sunday, serial
      integer :: epi_year_value

      serial = date_to_day_number(value)
      epi_year_value = epiyear_date(value)
      first_sunday = date_to_day_number(date_type(epi_year_value, 1, 4))
      first_sunday = first_sunday - int(wday_date(date_type(epi_year_value, 1, 4)) - 1, int64)
      component = int((serial - first_sunday) / 7_int64 + 1_int64)
   end function epiweek_date

   pure elemental integer function epiweek_datetime(value) result(component)
      !! Returns a date-time's epidemiological week number.
      type(datetime_type), intent(in) :: value !! Date-time to inspect.
      component = epiweek_date(value%date)
   end function epiweek_datetime

   pure elemental integer function epiyear_date(value) result(component)
      !! Returns the CDC epidemiological week-numbering year.
      type(date_type), intent(in) :: value !! Date to inspect.
      type(date_type) :: wednesday

      wednesday = date_from_day_number(date_to_day_number(value) + int(4 - wday_date(value), int64))
      component = wednesday%year
   end function epiyear_date

   pure elemental integer function epiyear_datetime(value) result(component)
      !! Returns a date-time's epidemiological week-numbering year.
      type(datetime_type), intent(in) :: value !! Date-time to inspect.
      component = epiyear_date(value%date)
   end function epiyear_datetime

   pure elemental integer function quarter_date(value, fiscal_start) result(component)
      !! Returns the quarter number relative to a fiscal starting month.
      type(date_type), intent(in) :: value !! Date to inspect.
      integer, intent(in), optional :: fiscal_start !! First fiscal month; defaults January.
      integer :: start

      start = 1
      if (present(fiscal_start)) start = fiscal_start
      component = modulo(value%month - start, 12) / 3 + 1
   end function quarter_date

   pure elemental integer function quarter_datetime(value, fiscal_start) result(component)
      !! Returns a date-time's fiscal quarter.
      type(datetime_type), intent(in) :: value !! Date-time to inspect.
      integer, intent(in), optional :: fiscal_start !! First fiscal month.
      component = quarter_date(value%date, fiscal_start)
   end function quarter_datetime

   pure elemental integer function semester_date(value, fiscal_start) result(component)
      !! Returns the half-year number relative to a fiscal starting month.
      type(date_type), intent(in) :: value !! Date to inspect.
      integer, intent(in), optional :: fiscal_start !! First fiscal month.
      integer :: start

      start = 1
      if (present(fiscal_start)) start = fiscal_start
      component = modulo(value%month - start, 12) / 6 + 1
   end function semester_date

   pure elemental integer function semester_datetime(value, fiscal_start) result(component)
      !! Returns a date-time's fiscal semester.
      type(datetime_type), intent(in) :: value !! Date-time to inspect.
      integer, intent(in), optional :: fiscal_start !! First fiscal month.
      component = semester_date(value%date, fiscal_start)
   end function semester_datetime

   pure elemental integer function days_in_month_date(value) result(number_of_days)
      !! Returns the number of days in a date's month.
      type(date_type), intent(in) :: value !! Date to inspect.
      number_of_days = days_in_month_components(value%year, value%month)
   end function days_in_month_date

   pure elemental integer function days_in_month_datetime(value) result(number_of_days)
      !! Returns the number of days in a date-time's month.
      type(datetime_type), intent(in) :: value !! Date-time to inspect.
      number_of_days = days_in_month_date(value%date)
   end function days_in_month_datetime

   pure elemental integer function days_in_month_components(year, month) result(number_of_days)
      !! Returns a Gregorian month length from numeric components.
      integer, intent(in) :: year, month !! Gregorian year and month.
      integer, parameter :: lengths(12) = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

      if (month < 1 .or. month > 12) then
         number_of_days = 0
      else
         number_of_days = lengths(month)
         if (month == 2 .and. leap_year(year)) number_of_days = 29
      end if
   end function days_in_month_components

   pure elemental logical function leap_year(year) result(is_leap)
      !! Reports whether a Gregorian year contains February 29.
      integer, intent(in) :: year !! Year to inspect.
      is_leap = mod(year, 4) == 0 .and. (mod(year, 100) /= 0 .or. mod(year, 400) == 0)
   end function leap_year

   pure elemental real(dp) function decimal_date_date(value) result(decimal)
      !! Converts a date to a fractional Gregorian year.
      type(date_type), intent(in) :: value !! Date to convert.
      integer :: length_of_year

      length_of_year = merge(366, 365, leap_year(value%year))
      decimal = real(value%year, dp) + real(yday_date(value) - 1, dp) / real(length_of_year, dp)
   end function decimal_date_date

   pure elemental real(dp) function decimal_date_datetime(value) result(decimal)
      !! Converts a date-time to a fractional Gregorian year.
      type(datetime_type), intent(in) :: value !! Date-time to convert.
      integer :: length_of_year

      length_of_year = merge(366, 365, leap_year(value%date%year))
      decimal = decimal_date_date(value%date)
      decimal = decimal + (real(value%hour * 3600 + value%minute * 60, dp) + value%second) / &
         (86400.0_dp * real(length_of_year, dp))
   end function decimal_date_datetime

   pure elemental function date_decimal(decimal) result(value)
      !! Converts a fractional Gregorian year to its containing civil date.
      real(dp), intent(in) :: decimal !! Fractional year.
      type(date_type) :: value
      integer :: day_index, year_value

      year_value = int(floor(decimal))
      day_index = int(floor((decimal - real(year_value, dp)) * &
         real(merge(366, 365, leap_year(year_value)), dp)))
      value = date_from_day_number(date_to_day_number(date_type(year_value, 1, 1)) + int(day_index, int64))
   end function date_decimal

   pure elemental function date_add_days(value, number_of_days) result(output)
      !! Adds whole days to a date.
      type(date_type), intent(in) :: value !! Starting date.
      integer, intent(in) :: number_of_days !! Signed day count.
      type(date_type) :: output

      output = date_from_day_number(date_to_day_number(value) + int(number_of_days, int64))
   end function date_add_days

   pure elemental function days_add_date(number_of_days, value) result(output)
      !! Adds whole days with the numeric operand first.
      integer, intent(in) :: number_of_days !! Signed day count.
      type(date_type), intent(in) :: value   !! Starting date.
      type(date_type) :: output
      output = date_add_days(value, number_of_days)
   end function days_add_date

   pure elemental function date_subtract_days(value, number_of_days) result(output)
      !! Subtracts whole days from a date.
      type(date_type), intent(in) :: value !! Starting date.
      integer, intent(in) :: number_of_days !! Signed day count.
      type(date_type) :: output
      output = date_add_days(value, -number_of_days)
   end function date_subtract_days

   pure elemental integer(int64) function date_difference(left, right) result(number_of_days)
      !! Returns the signed whole-day difference between dates.
      type(date_type), intent(in) :: left, right !! Dates to subtract.
      number_of_days = date_to_day_number(left) - date_to_day_number(right)
   end function date_difference

   pure elemental logical function date_equal(left, right) result(equal)
      !! Tests civil date equality.
      type(date_type), intent(in) :: left, right !! Dates to compare.
      equal = left%year == right%year .and. left%month == right%month .and. left%day == right%day
   end function date_equal

   pure elemental logical function date_not_equal(left, right) result(not_equal)
      !! Tests civil date inequality.
      type(date_type), intent(in) :: left, right !! Dates to compare.
      not_equal = .not. date_equal(left, right)
   end function date_not_equal

   pure elemental logical function date_less(left, right) result(less)
      !! Tests chronological date ordering.
      type(date_type), intent(in) :: left, right !! Dates to compare.
      less = date_to_day_number(left) < date_to_day_number(right)
   end function date_less

   pure elemental logical function date_less_equal(left, right) result(less_equal)
      !! Tests chronological date ordering with equality.
      type(date_type), intent(in) :: left, right !! Dates to compare.
      less_equal = date_to_day_number(left) <= date_to_day_number(right)
   end function date_less_equal

   pure elemental logical function date_greater(left, right) result(greater)
      !! Tests reverse chronological date ordering.
      type(date_type), intent(in) :: left, right !! Dates to compare.
      greater = date_to_day_number(left) > date_to_day_number(right)
   end function date_greater

   pure elemental logical function date_greater_equal(left, right) result(greater_equal)
      !! Tests reverse chronological date ordering with equality.
      type(date_type), intent(in) :: left, right !! Dates to compare.
      greater_equal = date_to_day_number(left) >= date_to_day_number(right)
   end function date_greater_equal

   pure elemental logical function datetime_equal(left, right) result(equal)
      !! Tests date-time instant equality across fixed offsets.
      type(datetime_type), intent(in) :: left, right !! Date-times to compare.
      equal = abs(datetime_to_epoch(left) - datetime_to_epoch(right)) <= 1.0e-9_dp
   end function datetime_equal

   pure elemental logical function datetime_not_equal(left, right) result(not_equal)
      !! Tests date-time instant inequality.
      type(datetime_type), intent(in) :: left, right !! Date-times to compare.
      not_equal = .not. datetime_equal(left, right)
   end function datetime_not_equal

   pure elemental logical function datetime_less(left, right) result(less)
      !! Tests chronological date-time ordering.
      type(datetime_type), intent(in) :: left, right !! Date-times to compare.
      less = datetime_to_epoch(left) < datetime_to_epoch(right)
   end function datetime_less

   pure elemental logical function datetime_less_equal(left, right) result(less_equal)
      !! Tests chronological date-time ordering with equality.
      type(datetime_type), intent(in) :: left, right !! Date-times to compare.
      less_equal = datetime_to_epoch(left) <= datetime_to_epoch(right)
   end function datetime_less_equal

   pure elemental logical function datetime_greater(left, right) result(greater)
      !! Tests reverse chronological date-time ordering.
      type(datetime_type), intent(in) :: left, right !! Date-times to compare.
      greater = datetime_to_epoch(left) > datetime_to_epoch(right)
   end function datetime_greater

   pure elemental logical function datetime_greater_equal(left, right) result(greater_equal)
      !! Tests reverse chronological date-time ordering with equality.
      type(datetime_type), intent(in) :: left, right !! Date-times to compare.
      greater_equal = datetime_to_epoch(left) >= datetime_to_epoch(right)
   end function datetime_greater_equal

   pure elemental integer(int64) function floor_divide(numerator, denominator) result(quotient)
      !! Computes mathematical floor division for signed integers.
      integer(int64), intent(in) :: numerator, denominator !! Dividend and positive divisor.

      quotient = numerator / denominator
      if (mod(numerator, denominator) < 0_int64) quotient = quotient - 1_int64
   end function floor_divide

end module lubridate_calendar
