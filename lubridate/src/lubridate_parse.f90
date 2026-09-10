! SPDX-License-Identifier: MIT
! SPDX-FileComment: Numeric parsers for the Fortran lubridate translation.
module lubridate_parse
   !! Parses common numeric date and date-time representations.
   use lubridate_calendar, only: make_date, make_datetime
   use lubridate_types, only: date_type, datetime_type, dp, period_type
   implicit none
   private

   public :: dmy, dmy_h, dmy_hm, dmy_hms, dym, hm, hms, mdy, mdy_h, mdy_hm, mdy_hms, ms, my, myd
   public :: parse_date_time, ydm, ydm_h, ydm_hm, ydm_hms, ym, ymd, ymd_h, ymd_hm, ymd_hms, yq

contains

   pure elemental function ymd(text) result(value)
      !! Parses a year-month-day string with numeric components.
      character(len=*), intent(in) :: text !! Text to parse.
      type(date_type) :: value
      value = parse_date_order(text, 'ymd')
   end function ymd

   pure elemental function ydm(text) result(value)
      !! Parses a year-day-month string with numeric components.
      character(len=*), intent(in) :: text !! Text to parse.
      type(date_type) :: value
      value = parse_date_order(text, 'ydm')
   end function ydm

   pure elemental function mdy(text) result(value)
      !! Parses a month-day-year string with numeric components.
      character(len=*), intent(in) :: text !! Text to parse.
      type(date_type) :: value
      value = parse_date_order(text, 'mdy')
   end function mdy

   pure elemental function myd(text) result(value)
      !! Parses a month-year-day string with numeric components.
      character(len=*), intent(in) :: text !! Text to parse.
      type(date_type) :: value
      value = parse_date_order(text, 'myd')
   end function myd

   pure elemental function dmy(text) result(value)
      !! Parses a day-month-year string with numeric components.
      character(len=*), intent(in) :: text !! Text to parse.
      type(date_type) :: value
      value = parse_date_order(text, 'dmy')
   end function dmy

   pure elemental function dym(text) result(value)
      !! Parses a day-year-month string with numeric components.
      character(len=*), intent(in) :: text !! Text to parse.
      type(date_type) :: value
      value = parse_date_order(text, 'dym')
   end function dym

   pure elemental function ym(text) result(value)
      !! Parses a year-month string, using the first day.
      character(len=*), intent(in) :: text !! Text to parse.
      type(date_type) :: value
      integer :: fields(8), count
      call numeric_fields(text, fields, count)
      if (count >= 2) value = make_date(fields(1), fields(2), 1)
   end function ym

   pure elemental function my(text) result(value)
      !! Parses a month-year string, using the first day.
      character(len=*), intent(in) :: text !! Text to parse.
      type(date_type) :: value
      integer :: fields(8), count
      call numeric_fields(text, fields, count)
      if (count >= 2) value = make_date(fields(2), fields(1), 1)
   end function my

   pure elemental function yq(text) result(value)
      !! Parses a year-quarter string as the quarter's first day.
      character(len=*), intent(in) :: text !! Text to parse.
      type(date_type) :: value
      integer :: fields(8), count
      call numeric_fields(text, fields, count)
      if (count >= 2 .and. fields(2) >= 1 .and. fields(2) <= 4) then
         value = make_date(fields(1), 3 * fields(2) - 2, 1)
      end if
   end function yq

   pure elemental function ymd_hms(text) result(value)
      !! Parses year, month, day, hour, minute, and second components.
      character(len=*), intent(in) :: text !! Text to parse.
      type(datetime_type) :: value
      value = parse_datetime_order(text, 'ymd', 3)
   end function ymd_hms

   pure elemental function ymd_hm(text) result(value)
      !! Parses year, month, day, hour, and minute components.
      character(len=*), intent(in) :: text !! Text to parse.
      type(datetime_type) :: value
      value = parse_datetime_order(text, 'ymd', 2)
   end function ymd_hm

   pure elemental function ymd_h(text) result(value)
      !! Parses year, month, day, and hour components.
      character(len=*), intent(in) :: text !! Text to parse.
      type(datetime_type) :: value
      value = parse_datetime_order(text, 'ymd', 1)
   end function ymd_h

   pure elemental function mdy_hms(text) result(value)
      !! Parses month-day-year followed by a time.
      character(len=*), intent(in) :: text !! Text to parse.
      type(datetime_type) :: value
      value = parse_datetime_order(text, 'mdy', 3)
   end function mdy_hms

   pure elemental function mdy_hm(text) result(value)
      !! Parses month-day-year followed by hour and minute.
      character(len=*), intent(in) :: text !! Text to parse.
      type(datetime_type) :: value
      value = parse_datetime_order(text, 'mdy', 2)
   end function mdy_hm

   pure elemental function mdy_h(text) result(value)
      !! Parses month-day-year followed by an hour.
      character(len=*), intent(in) :: text !! Text to parse.
      type(datetime_type) :: value
      value = parse_datetime_order(text, 'mdy', 1)
   end function mdy_h

   pure elemental function dmy_hms(text) result(value)
      !! Parses day-month-year followed by a time.
      character(len=*), intent(in) :: text !! Text to parse.
      type(datetime_type) :: value
      value = parse_datetime_order(text, 'dmy', 3)
   end function dmy_hms

   pure elemental function dmy_hm(text) result(value)
      !! Parses day-month-year followed by hour and minute.
      character(len=*), intent(in) :: text !! Text to parse.
      type(datetime_type) :: value
      value = parse_datetime_order(text, 'dmy', 2)
   end function dmy_hm

   pure elemental function dmy_h(text) result(value)
      !! Parses day-month-year followed by an hour.
      character(len=*), intent(in) :: text !! Text to parse.
      type(datetime_type) :: value
      value = parse_datetime_order(text, 'dmy', 1)
   end function dmy_h

   pure elemental function ydm_hms(text) result(value)
      !! Parses year-day-month followed by a time.
      character(len=*), intent(in) :: text !! Text to parse.
      type(datetime_type) :: value
      value = parse_datetime_order(text, 'ydm', 3)
   end function ydm_hms

   pure elemental function ydm_hm(text) result(value)
      !! Parses year-day-month followed by hour and minute.
      character(len=*), intent(in) :: text !! Text to parse.
      type(datetime_type) :: value
      value = parse_datetime_order(text, 'ydm', 2)
   end function ydm_hm

   pure elemental function ydm_h(text) result(value)
      !! Parses year-day-month followed by an hour.
      character(len=*), intent(in) :: text !! Text to parse.
      type(datetime_type) :: value
      value = parse_datetime_order(text, 'ydm', 1)
   end function ydm_h

   pure elemental function parse_date_time(text, order) result(value)
      !! Parses a date-time using a three-letter date order such as `ymd`.
      character(len=*), intent(in) :: text  !! Text to parse.
      character(len=*), intent(in) :: order !! Date component order.
      type(datetime_type) :: value
      value = parse_datetime_order(text, order, 3)
   end function parse_date_time

   pure elemental function hms(text) result(value)
      !! Parses hour-minute-second text as a calendar period.
      character(len=*), intent(in) :: text !! Text to parse.
      type(period_type) :: value
      integer :: fields(8), count
      call numeric_fields(text, fields, count)
      if (count >= 1) value%hours = fields(1)
      if (count >= 2) value%minutes = fields(2)
      if (count >= 3) value%seconds = real(fields(3), dp)
   end function hms

   pure elemental function hm(text) result(value)
      !! Parses hour-minute text as a calendar period.
      character(len=*), intent(in) :: text !! Text to parse.
      type(period_type) :: value
      integer :: fields(8), count
      call numeric_fields(text, fields, count)
      if (count >= 1) value%hours = fields(1)
      if (count >= 2) value%minutes = fields(2)
   end function hm

   pure elemental function ms(text) result(value)
      !! Parses minute-second text as a calendar period.
      character(len=*), intent(in) :: text !! Text to parse.
      type(period_type) :: value
      integer :: fields(8), count
      call numeric_fields(text, fields, count)
      if (count >= 1) value%minutes = fields(1)
      if (count >= 2) value%seconds = real(fields(2), dp)
   end function ms

   pure elemental function parse_date_order(text, order) result(value)
      !! Parses three date fields according to an explicit order.
      character(len=*), intent(in) :: text  !! Text to parse.
      character(len=*), intent(in) :: order !! Permutation of `ymd`.
      type(date_type) :: value
      integer :: fields(8), count, components(3), i
      call numeric_fields(text, fields, count)
      if (count < 3) call compact_fields(text, order, fields, count)
      if (count < 3 .or. len_trim(order) < 3) return
      components = 0
      do i = 1, 3
         select case (order(i:i))
         case ('y', 'Y')
            components(1) = fields(i)
         case ('m', 'M')
            components(2) = fields(i)
         case ('d', 'D')
            components(3) = fields(i)
         end select
      end do
      value = make_date(components(1), components(2), components(3))
   end function parse_date_order

   pure elemental function parse_datetime_order(text, order, time_fields) result(value)
      !! Parses three date fields followed by up to three time fields.
      character(len=*), intent(in) :: text  !! Text to parse.
      character(len=*), intent(in) :: order !! Permutation of `ymd`.
      integer, intent(in) :: time_fields    !! Number of expected time fields.
      type(datetime_type) :: value
      integer :: fields(8), count, components(3), i, hh, mm
      real(dp) :: ss
      call numeric_fields(text, fields, count)
      if (count < 6) call compact_fields(text, order, fields, count)
      if (count < 3 .or. len_trim(order) < 3) return
      components = 0
      do i = 1, 3
         select case (order(i:i))
         case ('y', 'Y')
            components(1) = fields(i)
         case ('m', 'M')
            components(2) = fields(i)
         case ('d', 'D')
            components(3) = fields(i)
         end select
      end do
      hh = 0
      mm = 0
      ss = 0.0_dp
      if (time_fields >= 1 .and. count >= 4) hh = fields(4)
      if (time_fields >= 2 .and. count >= 5) mm = fields(5)
      if (time_fields >= 3 .and. count >= 6) then
         ss = real(fields(6), dp)
         if (count >= 7) ss = ss + decimal_fraction(text)
      end if
      value = make_datetime(components(1), components(2), components(3), hh, mm, ss)
   end function parse_datetime_order

   pure subroutine numeric_fields(text, fields, count)
      !! Extracts unsigned integer fields separated by non-digits.
      character(len=*), intent(in) :: text !! Text containing numeric fields.
      integer, intent(out) :: fields(:)    !! Extracted values.
      integer, intent(out) :: count        !! Number of extracted fields.
      integer :: i, digit
      logical :: in_field
      fields = 0
      count = 0
      in_field = .false.
      do i = 1, len_trim(text)
         digit = iachar(text(i:i)) - iachar('0')
         if (digit >= 0 .and. digit <= 9) then
            if (.not. in_field) then
               if (count == size(fields)) return
               count = count + 1
               in_field = .true.
            end if
            fields(count) = 10 * fields(count) + digit
         else
            in_field = .false.
         end if
      end do
   end subroutine numeric_fields

   pure subroutine compact_fields(text, order, fields, count)
      !! Splits compact four-digit-year date and time representations.
      character(len=*), intent(in) :: text  !! Text containing compact digits.
      character(len=*), intent(in) :: order !! Permutation of `ymd`.
      integer, intent(out) :: fields(:)     !! Extracted ordered fields.
      integer, intent(out) :: count         !! Number of extracted fields.
      character(len=32) :: digits
      integer :: date_widths(3), digit_count, i, position, width

      digits = ''
      digit_count = 0
      do i = 1, len_trim(text)
         if (text(i:i) >= '0' .and. text(i:i) <= '9') then
            if (digit_count == len(digits)) exit
            digit_count = digit_count + 1
            digits(digit_count:digit_count) = text(i:i)
         end if
      end do
      fields = 0
      count = 0
      if (digit_count < 8 .or. len_trim(order) < 3) return
      do i = 1, 3
         date_widths(i) = merge(4, 2, order(i:i) == 'y' .or. order(i:i) == 'Y')
      end do
      position = 1
      do i = 1, 3
         width = date_widths(i)
         read (digits(position:position + width - 1), '(i4)') fields(i)
         position = position + width
      end do
      count = 3
      do while (digit_count - position + 1 >= 2 .and. count < min(6, size(fields)))
         count = count + 1
         read (digits(position:position + 1), '(i2)') fields(count)
         position = position + 2
      end do
   end subroutine compact_fields

   pure elemental real(dp) function decimal_fraction(text) result(fraction)
      !! Extracts digits following the final decimal point as a fraction.
      character(len=*), intent(in) :: text !! Text that may contain fractional seconds.
      integer :: decimal_position, digit, i, scale
      fraction = 0.0_dp
      decimal_position = 0
      do i = 1, len_trim(text)
         if (text(i:i) == '.') decimal_position = i
      end do
      if (decimal_position == 0) return
      scale = 1
      do i = decimal_position + 1, len_trim(text)
         digit = iachar(text(i:i)) - iachar('0')
         if (digit < 0 .or. digit > 9) exit
         scale = 10 * scale
         fraction = fraction + real(digit, dp) / real(scale, dp)
      end do
   end function decimal_fraction

end module lubridate_parse
