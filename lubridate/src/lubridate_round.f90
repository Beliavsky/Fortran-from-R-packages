! SPDX-License-Identifier: MIT
! SPDX-FileComment: Rounding and zone utilities for the Fortran lubridate translation.
module lubridate_round
   !! Implements calendar rounding, fixed-offset conversion, and clock queries.
   use, intrinsic :: iso_fortran_env, only: int64
   use lubridate_calendar, only: date_from_day_number, date_to_day_number
   use lubridate_calendar, only: datetime_from_epoch, datetime_to_epoch, make_date, make_datetime
   use lubridate_calendar, only: wday, operator(/=)
   use lubridate_spans, only: days, months, years, operator(+)
   use lubridate_types, only: cyclic_encoding_type, date_type, datetime_type, dp
   implicit none
   private

   public :: am, ceiling_date, cyclic_encoding, floor_date, force_tz, now, pm
   public :: round_date, today, with_tz

   interface floor_date
      module procedure floor_date_value, floor_datetime_value
   end interface floor_date

   interface ceiling_date
      module procedure ceiling_date_value, ceiling_datetime_value
   end interface ceiling_date

   interface round_date
      module procedure round_date_value, round_datetime_value
   end interface round_date

contains

   pure elemental function floor_date_value(value, unit, week_start) result(output)
      !! Rounds a date down to a calendar boundary.
      type(date_type), intent(in) :: value        !! Date to round.
      character(len=*), intent(in) :: unit        !! Calendar unit.
      integer, intent(in), optional :: week_start !! First weekday, Sunday equals one.
      type(date_type) :: output
      integer :: start
      start = 1
      if (present(week_start)) start = week_start
      select case (trim(adjustl(unit)))
      case ('day', 'days')
         output = value
      case ('week', 'weeks')
         output = date_from_day_number(date_to_day_number(value) - int(wday(value, start) - 1, int64))
      case ('month', 'months')
         output = make_date(value%year, value%month, 1)
      case ('quarter', 'quarters')
         output = make_date(value%year, 3 * ((value%month - 1) / 3) + 1, 1)
      case ('year', 'years')
         output = make_date(value%year, 1, 1)
      case default
         output = value
      end select
   end function floor_date_value

   pure elemental function ceiling_date_value(value, unit, week_start, change_on_boundary) result(output)
      !! Rounds a date up to a calendar boundary.
      type(date_type), intent(in) :: value        !! Date to round.
      character(len=*), intent(in) :: unit        !! Calendar unit.
      integer, intent(in), optional :: week_start !! First weekday, Sunday equals one.
      logical, intent(in), optional :: change_on_boundary !! Advance exact boundaries when true.
      type(date_type) :: output, lower
      logical :: advance
      integer :: amount
      lower = floor_date_value(value, unit, week_start)
      advance = value /= lower
      if (present(change_on_boundary)) advance = advance .or. change_on_boundary
      if (.not. advance) then
         output = lower
         return
      end if
      select case (trim(adjustl(unit)))
      case ('day', 'days')
         output = lower + days(1)
      case ('week', 'weeks')
         output = lower + days(7)
      case ('month', 'months')
         output = lower + months(1)
      case ('quarter', 'quarters')
         output = lower + months(3)
      case ('year', 'years')
         output = lower + years(1)
      case default
         amount = 0
         output = lower + days(amount)
      end select
   end function ceiling_date_value

   pure elemental function round_date_value(value, unit, week_start) result(output)
      !! Rounds a date to the nearest calendar boundary, breaking ties upward.
      type(date_type), intent(in) :: value        !! Date to round.
      character(len=*), intent(in) :: unit        !! Calendar unit.
      integer, intent(in), optional :: week_start !! First weekday, Sunday equals one.
      type(date_type) :: output, lower, upper
      lower = floor_date_value(value, unit, week_start)
      upper = ceiling_date_value(value, unit, week_start, change_on_boundary=.true.)
      if (date_to_day_number(value) - date_to_day_number(lower) < &
          date_to_day_number(upper) - date_to_day_number(value)) then
         output = lower
      else
         output = upper
      end if
   end function round_date_value

   pure elemental function floor_datetime_value(value, unit, week_start) result(output)
      !! Rounds a date-time down in its displayed fixed offset.
      type(datetime_type), intent(in) :: value    !! Date-time to round.
      character(len=*), intent(in) :: unit        !! Calendar or clock unit.
      integer, intent(in), optional :: week_start !! First weekday, Sunday equals one.
      type(datetime_type) :: output
      type(date_type) :: rounded_date
      real(dp) :: local_seconds, factor, rounded
      select case (trim(adjustl(unit)))
      case ('second', 'seconds', 'minute', 'minutes', 'hour', 'hours')
         factor = clock_unit_seconds(unit)
         local_seconds = datetime_to_epoch(value) + 60.0_dp * real(value%offset_minutes, dp)
         rounded = floor(local_seconds / factor) * factor
         output = datetime_from_epoch(rounded - 60.0_dp * real(value%offset_minutes, dp), value%offset_minutes)
      case default
         rounded_date = floor_date_value(value%date, unit, week_start)
         output = make_datetime(rounded_date%year, rounded_date%month, rounded_date%day, &
                                0, 0, 0.0_dp, value%offset_minutes)
      end select
   end function floor_datetime_value

   pure elemental function ceiling_datetime_value(value, unit, week_start, change_on_boundary) result(output)
      !! Rounds a date-time up in its displayed fixed offset.
      type(datetime_type), intent(in) :: value    !! Date-time to round.
      character(len=*), intent(in) :: unit        !! Calendar or clock unit.
      integer, intent(in), optional :: week_start !! First weekday, Sunday equals one.
      logical, intent(in), optional :: change_on_boundary !! Advance exact boundaries when true.
      type(datetime_type) :: output, lower
      type(date_type) :: next_date
      real(dp) :: factor
      logical :: advance
      lower = floor_datetime_value(value, unit, week_start)
      advance = abs(datetime_to_epoch(value) - datetime_to_epoch(lower)) > 1.0e-9_dp
      if (present(change_on_boundary)) advance = advance .or. change_on_boundary
      if (.not. advance) then
         output = lower
         return
      end if
      select case (trim(adjustl(unit)))
      case ('second', 'seconds', 'minute', 'minutes', 'hour', 'hours')
         factor = clock_unit_seconds(unit)
         output = datetime_from_epoch(datetime_to_epoch(lower) + factor, lower%offset_minutes)
      case default
         next_date = ceiling_date_value(value%date, unit, week_start, change_on_boundary=.true.)
         output = make_datetime(next_date%year, next_date%month, next_date%day, 0, 0, 0.0_dp, value%offset_minutes)
      end select
   end function ceiling_datetime_value

   pure elemental function round_datetime_value(value, unit, week_start) result(output)
      !! Rounds a date-time to the nearest boundary, breaking ties upward.
      type(datetime_type), intent(in) :: value    !! Date-time to round.
      character(len=*), intent(in) :: unit        !! Calendar or clock unit.
      integer, intent(in), optional :: week_start !! First weekday, Sunday equals one.
      type(datetime_type) :: output, lower, upper
      lower = floor_datetime_value(value, unit, week_start)
      upper = ceiling_datetime_value(value, unit, week_start, change_on_boundary=.true.)
      if (datetime_to_epoch(value) - datetime_to_epoch(lower) < &
          datetime_to_epoch(upper) - datetime_to_epoch(value)) then
         output = lower
      else
         output = upper
      end if
   end function round_datetime_value

   pure elemental function with_tz(value, offset_minutes) result(output)
      !! Displays the same instant using another fixed UTC offset.
      type(datetime_type), intent(in) :: value !! Date-time to convert.
      integer, intent(in) :: offset_minutes    !! New minutes east of UTC.
      type(datetime_type) :: output
      output = datetime_from_epoch(datetime_to_epoch(value), offset_minutes)
   end function with_tz

   pure elemental function force_tz(value, offset_minutes) result(output)
      !! Changes a fixed offset while preserving displayed clock components.
      type(datetime_type), intent(in) :: value !! Date-time to relabel.
      integer, intent(in) :: offset_minutes    !! Replacement minutes east of UTC.
      type(datetime_type) :: output
      output = value
      output%offset_minutes = offset_minutes
   end function force_tz

   pure elemental logical function am(value) result(is_am)
      !! Tests whether a displayed time is before noon.
      type(datetime_type), intent(in) :: value !! Date-time to inspect.
      is_am = value%hour < 12
   end function am

   pure elemental logical function pm(value) result(is_pm)
      !! Tests whether a displayed time is noon or later.
      type(datetime_type), intent(in) :: value !! Date-time to inspect.
      is_pm = value%hour >= 12
   end function pm

   pure elemental function cyclic_encoding(value, period, offset) result(encoded)
      !! Encodes a periodic scalar as sine and cosine coordinates.
      real(dp), intent(in) :: value           !! Value on the cycle.
      real(dp), intent(in) :: period          !! Cycle length.
      real(dp), intent(in), optional :: offset !! Cycle origin.
      type(cyclic_encoding_type) :: encoded
      real(dp) :: angle, origin_value
      origin_value = 0.0_dp
      if (present(offset)) origin_value = offset
      angle = 2.0_dp * acos(-1.0_dp) * (value - origin_value) / period
      encoded%sine = sin(angle)
      encoded%cosine = cos(angle)
   end function cyclic_encoding

   function now(offset_minutes) result(value)
      !! Returns the current system clock as a date-time.
      integer, intent(in), optional :: offset_minutes !! Requested fixed offset.
      type(datetime_type) :: value
      integer :: fields(8), offset
      call date_and_time(values=fields)
      offset = fields(4)
      value = make_datetime(fields(1), fields(2), fields(3), fields(5), fields(6), &
                            real(fields(7), dp) + real(fields(8), dp) / 1000.0_dp, offset)
      if (present(offset_minutes)) value = with_tz(value, offset_minutes)
   end function now

   function today(offset_minutes) result(value)
      !! Returns the current date in a requested fixed offset.
      integer, intent(in), optional :: offset_minutes !! Requested fixed offset.
      type(date_type) :: value
      type(datetime_type) :: current
      current = now(offset_minutes)
      value = current%date
   end function today

   pure elemental real(dp) function clock_unit_seconds(unit) result(factor)
      !! Returns seconds in a supported clock unit.
      character(len=*), intent(in) :: unit !! Unit name.
      select case (trim(adjustl(unit)))
      case ('minute', 'minutes')
         factor = 60.0_dp
      case ('hour', 'hours')
         factor = 3600.0_dp
      case default
         factor = 1.0_dp
      end select
   end function clock_unit_seconds

end module lubridate_round
