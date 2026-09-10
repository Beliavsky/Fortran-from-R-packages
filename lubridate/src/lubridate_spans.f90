! SPDX-License-Identifier: MIT
! SPDX-FileComment: Time-span arithmetic for the Fortran lubridate translation.
module lubridate_spans
   !! Implements exact durations, calendar periods, and time intervals.
   use lubridate_calendar, only: date_from_day_number, date_to_day_number
   use lubridate_calendar, only: datetime_from_epoch, datetime_to_epoch, days_in_month
   use lubridate_types, only: date_type, datetime_type, dp, duration_type, interval_type, period_type
   implicit none
   private

   public :: add_with_rollback, ddays, dhours, dmicroseconds, dmilliseconds, dminutes
   public :: dmonths, dnanoseconds, dpicoseconds, dseconds, duration
   public :: dweeks, dyears, int_aligns, int_end, int_flip, int_length, int_overlaps
   public :: int_shift, int_standardize, interval, period, period_to_seconds
   public :: rollback, rollforward, seconds_to_period, time_length, within
   public :: days, hours, microseconds, milliseconds, minutes, months, nanoseconds
   public :: picoseconds, seconds, weeks, years
   public :: operator(+), operator(-), operator(*), operator(/)

   interface operator(+)
      module procedure duration_add, date_add_period, period_add_date
      module procedure datetime_add_duration, duration_add_datetime
      module procedure datetime_add_period, period_add_datetime
   end interface operator(+)

   interface operator(-)
      module procedure duration_subtract, duration_negate, datetime_subtract_duration
      module procedure datetime_difference, date_subtract_period, datetime_subtract_period
   end interface operator(-)

   interface operator(*)
      module procedure duration_multiply_real, real_multiply_duration
   end interface operator(*)

   interface operator(/)
      module procedure duration_divide_real
   end interface operator(/)

   interface time_length
      module procedure time_length_duration, time_length_interval
   end interface time_length

contains

   pure elemental function duration(seconds) result(value)
      !! Constructs an exact duration measured in seconds.
      real(dp), intent(in) :: seconds !! Number of seconds.
      type(duration_type) :: value
      value%seconds = seconds
   end function duration

   pure elemental function dseconds(value) result(span)
      !! Constructs a duration from seconds.
      real(dp), intent(in) :: value !! Number of seconds.
      type(duration_type) :: span
      span = duration(value)
   end function dseconds

   pure elemental function dmilliseconds(value) result(span)
      !! Constructs a duration from milliseconds.
      real(dp), intent(in) :: value !! Number of milliseconds.
      type(duration_type) :: span
      span = duration(value * 1.0e-3_dp)
   end function dmilliseconds

   pure elemental function dmicroseconds(value) result(span)
      !! Constructs a duration from microseconds.
      real(dp), intent(in) :: value !! Number of microseconds.
      type(duration_type) :: span
      span = duration(value * 1.0e-6_dp)
   end function dmicroseconds

   pure elemental function dnanoseconds(value) result(span)
      !! Constructs a duration from nanoseconds.
      real(dp), intent(in) :: value !! Number of nanoseconds.
      type(duration_type) :: span
      span = duration(value * 1.0e-9_dp)
   end function dnanoseconds

   pure elemental function dpicoseconds(value) result(span)
      !! Constructs a duration from picoseconds.
      real(dp), intent(in) :: value !! Number of picoseconds.
      type(duration_type) :: span
      span = duration(value * 1.0e-12_dp)
   end function dpicoseconds

   pure elemental function dminutes(value) result(span)
      !! Constructs a duration from 60-second minutes.
      real(dp), intent(in) :: value !! Number of minutes.
      type(duration_type) :: span
      span = duration(60.0_dp * value)
   end function dminutes

   pure elemental function dhours(value) result(span)
      !! Constructs a duration from 3600-second hours.
      real(dp), intent(in) :: value !! Number of hours.
      type(duration_type) :: span
      span = duration(3600.0_dp * value)
   end function dhours

   pure elemental function ddays(value) result(span)
      !! Constructs a duration from 86400-second days.
      real(dp), intent(in) :: value !! Number of days.
      type(duration_type) :: span
      span = duration(86400.0_dp * value)
   end function ddays

   pure elemental function dweeks(value) result(span)
      !! Constructs a duration from seven-day weeks.
      real(dp), intent(in) :: value !! Number of weeks.
      type(duration_type) :: span
      span = ddays(7.0_dp * value)
   end function dweeks

   pure elemental function dmonths(value) result(span)
      !! Constructs a duration using lubridate's average Gregorian month.
      real(dp), intent(in) :: value !! Number of average months.
      type(duration_type) :: span
      span = duration(value * 31557600.0_dp / 12.0_dp)
   end function dmonths

   pure elemental function dyears(value) result(span)
      !! Constructs a duration using lubridate's 365.25-day year.
      real(dp), intent(in) :: value !! Number of average years.
      type(duration_type) :: span
      span = duration(value * 31557600.0_dp)
   end function dyears

   pure elemental function period(years, months, weeks, days, hours, minutes, seconds) result(value)
      !! Constructs a calendar period from optional components.
      integer, intent(in), optional :: years   !! Calendar years.
      integer, intent(in), optional :: months  !! Calendar months.
      integer, intent(in), optional :: weeks   !! Seven-day groups.
      integer, intent(in), optional :: days    !! Calendar days.
      integer, intent(in), optional :: hours   !! Hours.
      integer, intent(in), optional :: minutes !! Minutes.
      real(dp), intent(in), optional :: seconds !! Seconds.
      type(period_type) :: value
      if (present(years)) value%years = years
      if (present(months)) value%months = months
      if (present(weeks)) value%days = 7 * weeks
      if (present(days)) value%days = value%days + days
      if (present(hours)) value%hours = hours
      if (present(minutes)) value%minutes = minutes
      if (present(seconds)) value%seconds = seconds
   end function period

   pure elemental function seconds(value) result(span)
      !! Constructs a calendar period containing seconds.
      real(dp), intent(in) :: value !! Number of seconds.
      type(period_type) :: span
      span%seconds = value
   end function seconds

   pure elemental function milliseconds(value) result(span)
      !! Constructs a calendar period from milliseconds.
      real(dp), intent(in) :: value !! Number of milliseconds.
      type(period_type) :: span
      span%seconds = value * 1.0e-3_dp
   end function milliseconds

   pure elemental function microseconds(value) result(span)
      !! Constructs a calendar period from microseconds.
      real(dp), intent(in) :: value !! Number of microseconds.
      type(period_type) :: span
      span%seconds = value * 1.0e-6_dp
   end function microseconds

   pure elemental function nanoseconds(value) result(span)
      !! Constructs a calendar period from nanoseconds.
      real(dp), intent(in) :: value !! Number of nanoseconds.
      type(period_type) :: span
      span%seconds = value * 1.0e-9_dp
   end function nanoseconds

   pure elemental function picoseconds(value) result(span)
      !! Constructs a calendar period from picoseconds.
      real(dp), intent(in) :: value !! Number of picoseconds.
      type(period_type) :: span
      span%seconds = value * 1.0e-12_dp
   end function picoseconds

   pure elemental function minutes(value) result(span)
      !! Constructs a calendar period containing minutes.
      integer, intent(in) :: value !! Number of minutes.
      type(period_type) :: span
      span%minutes = value
   end function minutes

   pure elemental function hours(value) result(span)
      !! Constructs a calendar period containing hours.
      integer, intent(in) :: value !! Number of hours.
      type(period_type) :: span
      span%hours = value
   end function hours

   pure elemental function days(value) result(span)
      !! Constructs a calendar period containing days.
      integer, intent(in) :: value !! Number of days.
      type(period_type) :: span
      span%days = value
   end function days

   pure elemental function weeks(value) result(span)
      !! Constructs a calendar period containing seven-day weeks.
      integer, intent(in) :: value !! Number of weeks.
      type(period_type) :: span
      span%days = 7 * value
   end function weeks

   pure elemental function months(value) result(span)
      !! Constructs a calendar period containing months.
      integer, intent(in) :: value !! Number of months.
      type(period_type) :: span
      span%months = value
   end function months

   pure elemental function years(value) result(span)
      !! Constructs a calendar period containing years.
      integer, intent(in) :: value !! Number of years.
      type(period_type) :: span
      span%years = value
   end function years

   pure elemental real(dp) function period_to_seconds(value) result(total)
      !! Approximates a period in seconds using average Gregorian years and months.
      type(period_type), intent(in) :: value !! Period to convert.
      total = 31557600.0_dp * real(value%years, dp)
      total = total + 31557600.0_dp * real(value%months, dp) / 12.0_dp
      total = total + 86400.0_dp * real(value%days, dp)
      total = total + 3600.0_dp * real(value%hours, dp)
      total = total + 60.0_dp * real(value%minutes, dp) + value%seconds
   end function period_to_seconds

   pure elemental function seconds_to_period(value) result(span)
      !! Decomposes seconds into days, hours, minutes, and seconds.
      real(dp), intent(in) :: value !! Total seconds.
      type(period_type) :: span
      real(dp) :: remaining
      integer :: sign_value
      sign_value = merge(1, -1, value >= 0.0_dp)
      remaining = abs(value)
      span%days = sign_value * int(remaining / 86400.0_dp)
      remaining = modulo(remaining, 86400.0_dp)
      span%hours = sign_value * int(remaining / 3600.0_dp)
      remaining = modulo(remaining, 3600.0_dp)
      span%minutes = sign_value * int(remaining / 60.0_dp)
      span%seconds = real(sign_value, dp) * modulo(remaining, 60.0_dp)
   end function seconds_to_period

   pure elemental function interval(start, finish) result(value)
      !! Constructs a half-open time interval from two date-times.
      type(datetime_type), intent(in) :: start  !! Interval start.
      type(datetime_type), intent(in) :: finish !! Interval end.
      type(interval_type) :: value
      value%start = start
      value%finish = finish
   end function interval

   pure elemental function int_start(value) result(start)
      !! Returns an interval's start.
      type(interval_type), intent(in) :: value !! Interval to inspect.
      type(datetime_type) :: start
      start = value%start
   end function int_start

   pure elemental function int_end(value) result(finish)
      !! Returns an interval's end.
      type(interval_type), intent(in) :: value !! Interval to inspect.
      type(datetime_type) :: finish
      finish = value%finish
   end function int_end

   pure elemental real(dp) function int_length(value) result(length_seconds)
      !! Returns signed interval length in seconds.
      type(interval_type), intent(in) :: value !! Interval to measure.
      length_seconds = datetime_to_epoch(value%finish) - datetime_to_epoch(value%start)
   end function int_length

   pure elemental function int_flip(value) result(flipped)
      !! Reverses an interval.
      type(interval_type), intent(in) :: value !! Interval to reverse.
      type(interval_type) :: flipped
      flipped = interval(value%finish, value%start)
   end function int_flip

   pure elemental function int_standardize(value) result(standardized)
      !! Orders an interval so that its length is nonnegative.
      type(interval_type), intent(in) :: value !! Interval to standardize.
      type(interval_type) :: standardized
      if (int_length(value) < 0.0_dp) then
         standardized = int_flip(value)
      else
         standardized = value
      end if
   end function int_standardize

   pure elemental function int_shift(value, span) result(shifted)
      !! Shifts both endpoints by an exact duration.
      type(interval_type), intent(in) :: value !! Interval to shift.
      type(duration_type), intent(in) :: span   !! Exact shift.
      type(interval_type) :: shifted
      shifted = interval(value%start + span, value%finish + span)
   end function int_shift

   pure elemental logical function int_overlaps(left, right) result(overlaps)
      !! Tests whether two standardized closed intervals overlap.
      type(interval_type), intent(in) :: left  !! First interval.
      type(interval_type), intent(in) :: right !! Second interval.
      type(interval_type) :: a, b
      a = int_standardize(left)
      b = int_standardize(right)
      overlaps = datetime_to_epoch(a%start) <= datetime_to_epoch(b%finish) .and. &
                 datetime_to_epoch(b%start) <= datetime_to_epoch(a%finish)
   end function int_overlaps

   pure elemental logical function int_aligns(left, right) result(aligns)
      !! Tests whether two intervals share a start or end instant.
      type(interval_type), intent(in) :: left  !! First interval.
      type(interval_type), intent(in) :: right !! Second interval.
      aligns = datetime_to_epoch(left%start) == datetime_to_epoch(right%start) .or. &
               datetime_to_epoch(left%finish) == datetime_to_epoch(right%finish)
   end function int_aligns

   pure elemental logical function within(moment, value) result(is_within)
      !! Tests whether a date-time lies in a standardized closed interval.
      type(datetime_type), intent(in) :: moment !! Instant to test.
      type(interval_type), intent(in) :: value   !! Bounding interval.
      type(interval_type) :: standardized
      real(dp) :: instant
      standardized = int_standardize(value)
      instant = datetime_to_epoch(moment)
      is_within = instant >= datetime_to_epoch(standardized%start) .and. &
                  instant <= datetime_to_epoch(standardized%finish)
   end function within

   pure elemental real(dp) function time_length_duration(value, unit) result(length)
      !! Expresses a duration in a requested fixed unit.
      type(duration_type), intent(in) :: value !! Duration to measure.
      character(len=*), intent(in) :: unit     !! Unit name.
      length = value%seconds / unit_seconds(unit)
   end function time_length_duration

   pure elemental real(dp) function time_length_interval(value, unit) result(length)
      !! Expresses an interval length in a requested fixed unit.
      type(interval_type), intent(in) :: value !! Interval to measure.
      character(len=*), intent(in) :: unit     !! Unit name.
      length = int_length(value) / unit_seconds(unit)
   end function time_length_interval

   pure elemental function add_with_rollback(value, span, roll_to_first) result(output)
      !! Adds a period, clipping an invalid month day to a month boundary.
      type(date_type), intent(in) :: value      !! Date to adjust.
      type(period_type), intent(in) :: span     !! Calendar period to add.
      logical, intent(in), optional :: roll_to_first !! Roll overflow to next month when true.
      type(date_type) :: output
      integer :: total_month, target_year, target_month, target_day, month_days
      logical :: first
      first = .false.
      if (present(roll_to_first)) first = roll_to_first
      total_month = 12 * (value%year + span%years) + value%month - 1 + span%months
      target_year = floor_divide_int(total_month, 12)
      target_month = modulo(total_month, 12) + 1
      month_days = days_in_month(target_year, target_month)
      target_day = min(value%day, month_days)
      if (first .and. value%day > month_days) then
         total_month = total_month + 1
         target_year = floor_divide_int(total_month, 12)
         target_month = modulo(total_month, 12) + 1
         target_day = 1
      end if
      output = date_type(target_year, target_month, target_day)
      output = date_from_day_number(date_to_day_number(output) + int(span%days, kind(date_to_day_number(output))))
   end function add_with_rollback

   pure elemental function rollback(value) result(output)
      !! Returns the final day of the preceding month.
      type(date_type), intent(in) :: value !! Date whose month is rolled back.
      type(date_type) :: output
      output = date_from_day_number(date_to_day_number(date_type(value%year, value%month, 1)) - 1)
   end function rollback

   pure elemental function rollforward(value) result(output)
      !! Returns the first day of the following month.
      type(date_type), intent(in) :: value !! Date whose month is rolled forward.
      type(date_type) :: output
      output = add_with_rollback(date_type(value%year, value%month, 1), months(1))
   end function rollforward

   pure elemental function date_add_period(value, span) result(output)
      !! Adds a calendar period to a date with month rollback.
      type(date_type), intent(in) :: value  !! Date operand.
      type(period_type), intent(in) :: span !! Period operand.
      type(date_type) :: output
      output = add_with_rollback(value, span)
   end function date_add_period

   pure elemental function period_add_date(span, value) result(output)
      !! Adds a calendar period to a date.
      type(period_type), intent(in) :: span !! Period operand.
      type(date_type), intent(in) :: value  !! Date operand.
      type(date_type) :: output
      output = value + span
   end function period_add_date

   pure elemental function date_subtract_period(value, span) result(output)
      !! Subtracts a calendar period from a date.
      type(date_type), intent(in) :: value  !! Date operand.
      type(period_type), intent(in) :: span !! Period operand.
      type(date_type) :: output
      type(period_type) :: negative
      negative = negate_period(span)
      output = value + negative
   end function date_subtract_period

   pure elemental function datetime_add_duration(value, span) result(output)
      !! Adds an exact duration to a date-time instant.
      type(datetime_type), intent(in) :: value !! Date-time operand.
      type(duration_type), intent(in) :: span  !! Exact duration.
      type(datetime_type) :: output
      output = datetime_from_epoch(datetime_to_epoch(value) + span%seconds, value%offset_minutes)
   end function datetime_add_duration

   pure elemental function duration_add_datetime(span, value) result(output)
      !! Adds an exact duration to a date-time instant.
      type(duration_type), intent(in) :: span  !! Exact duration.
      type(datetime_type), intent(in) :: value !! Date-time operand.
      type(datetime_type) :: output
      output = value + span
   end function duration_add_datetime

   pure elemental function datetime_subtract_duration(value, span) result(output)
      !! Subtracts an exact duration from a date-time instant.
      type(datetime_type), intent(in) :: value !! Date-time operand.
      type(duration_type), intent(in) :: span  !! Exact duration.
      type(datetime_type) :: output
      output = datetime_from_epoch(datetime_to_epoch(value) - span%seconds, value%offset_minutes)
   end function datetime_subtract_duration

   pure elemental function datetime_difference(left, right) result(span)
      !! Returns the exact duration between two date-times.
      type(datetime_type), intent(in) :: left  !! Left operand.
      type(datetime_type), intent(in) :: right !! Right operand.
      type(duration_type) :: span
      span = duration(datetime_to_epoch(left) - datetime_to_epoch(right))
   end function datetime_difference

   pure elemental function datetime_add_period(value, span) result(output)
      !! Adds calendar and clock components to a date-time.
      type(datetime_type), intent(in) :: value !! Date-time operand.
      type(period_type), intent(in) :: span    !! Period operand.
      type(datetime_type) :: output
      type(date_type) :: adjusted_date
      real(dp) :: clock_seconds
      adjusted_date = value%date + period(years=span%years, months=span%months, days=span%days)
      output = datetime_type(adjusted_date, value%hour, value%minute, value%second, value%offset_minutes)
      clock_seconds = 3600.0_dp * real(span%hours, dp) + 60.0_dp * real(span%minutes, dp) + span%seconds
      output = output + duration(clock_seconds)
   end function datetime_add_period

   pure elemental function period_add_datetime(span, value) result(output)
      !! Adds calendar and clock components to a date-time.
      type(period_type), intent(in) :: span    !! Period operand.
      type(datetime_type), intent(in) :: value !! Date-time operand.
      type(datetime_type) :: output
      output = value + span
   end function period_add_datetime

   pure elemental function datetime_subtract_period(value, span) result(output)
      !! Subtracts calendar and clock components from a date-time.
      type(datetime_type), intent(in) :: value !! Date-time operand.
      type(period_type), intent(in) :: span    !! Period operand.
      type(datetime_type) :: output
      output = value + negate_period(span)
   end function datetime_subtract_period

   pure elemental function duration_add(left, right) result(output)
      !! Adds two exact durations.
      type(duration_type), intent(in) :: left  !! Left operand.
      type(duration_type), intent(in) :: right !! Right operand.
      type(duration_type) :: output
      output = duration(left%seconds + right%seconds)
   end function duration_add

   pure elemental function duration_subtract(left, right) result(output)
      !! Subtracts two exact durations.
      type(duration_type), intent(in) :: left  !! Left operand.
      type(duration_type), intent(in) :: right !! Right operand.
      type(duration_type) :: output
      output = duration(left%seconds - right%seconds)
   end function duration_subtract

   pure elemental function duration_negate(value) result(output)
      !! Negates an exact duration.
      type(duration_type), intent(in) :: value !! Duration to negate.
      type(duration_type) :: output
      output = duration(-value%seconds)
   end function duration_negate

   pure elemental function duration_multiply_real(value, factor) result(output)
      !! Multiplies an exact duration by a scalar.
      type(duration_type), intent(in) :: value !! Duration operand.
      real(dp), intent(in) :: factor           !! Scale factor.
      type(duration_type) :: output
      output = duration(value%seconds * factor)
   end function duration_multiply_real

   pure elemental function real_multiply_duration(factor, value) result(output)
      !! Multiplies an exact duration by a scalar.
      real(dp), intent(in) :: factor           !! Scale factor.
      type(duration_type), intent(in) :: value !! Duration operand.
      type(duration_type) :: output
      output = value * factor
   end function real_multiply_duration

   pure elemental function duration_divide_real(value, divisor) result(output)
      !! Divides an exact duration by a scalar.
      type(duration_type), intent(in) :: value !! Duration operand.
      real(dp), intent(in) :: divisor          !! Divisor.
      type(duration_type) :: output
      output = duration(value%seconds / divisor)
   end function duration_divide_real

   pure elemental function negate_period(value) result(output)
      !! Negates every component of a calendar period.
      type(period_type), intent(in) :: value !! Period to negate.
      type(period_type) :: output
      output = period_type(-value%years, -value%months, -value%days, -value%hours, &
                           -value%minutes, -value%seconds)
   end function negate_period

   pure elemental real(dp) function unit_seconds(unit) result(factor)
      !! Returns the fixed number of seconds in a supported unit.
      character(len=*), intent(in) :: unit !! Unit name.
      select case (trim(adjustl(unit)))
      case ('second', 'seconds')
         factor = 1.0_dp
      case ('minute', 'minutes')
         factor = 60.0_dp
      case ('hour', 'hours')
         factor = 3600.0_dp
      case ('day', 'days')
         factor = 86400.0_dp
      case ('week', 'weeks')
         factor = 604800.0_dp
      case ('month', 'months')
         factor = 31557600.0_dp / 12.0_dp
      case ('year', 'years')
         factor = 31557600.0_dp
      case default
         factor = 1.0_dp
      end select
   end function unit_seconds

   pure elemental integer function floor_divide_int(numerator, denominator) result(quotient)
      !! Computes mathematical floor division for integers.
      integer, intent(in) :: numerator   !! Dividend.
      integer, intent(in) :: denominator !! Positive divisor.
      quotient = numerator / denominator
      if (numerator < 0 .and. modulo(numerator, denominator) /= 0) quotient = quotient - 1
   end function floor_divide_int

end module lubridate_spans
