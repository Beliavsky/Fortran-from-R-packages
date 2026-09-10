! SPDX-License-Identifier: MIT
! SPDX-FileComment: Date-time types for the Fortran lubridate translation.
module lubridate_types
   !! Defines civil dates, fixed-offset date-times, and time-span values.
   use, intrinsic :: iso_fortran_env, only: real64
   implicit none
   private

   integer, parameter, public :: dp = real64

   type, public :: date_type
      !! Stores one proleptic-Gregorian civil date.
      integer :: year = 0
      integer :: month = 0
      integer :: day = 0
   contains
      procedure :: valid => date_valid
      procedure :: iso8601 => date_iso8601
   end type date_type

   type, public :: datetime_type
      !! Stores civil clock components with a fixed UTC offset in minutes.
      type(date_type) :: date
      integer :: hour = 0
      integer :: minute = 0
      real(dp) :: second = 0.0_dp
      integer :: offset_minutes = 0
   contains
      procedure :: valid => datetime_valid
      procedure :: iso8601 => datetime_iso8601
   end type datetime_type

   type, public :: duration_type
      !! Stores an exact physical duration in seconds.
      real(dp) :: seconds = 0.0_dp
   end type duration_type

   type, public :: period_type
      !! Stores a civil-time period whose components retain calendar meaning.
      integer :: years = 0
      integer :: months = 0
      integer :: days = 0
      integer :: hours = 0
      integer :: minutes = 0
      real(dp) :: seconds = 0.0_dp
   end type period_type

   type, public :: interval_type
      !! Stores a directed interval between two date-time instants.
      type(datetime_type) :: start
      type(datetime_type) :: finish
   end type interval_type

   type, public :: cyclic_encoding_type
      !! Stores sine and cosine coordinates for a cyclic component.
      real(dp) :: sine = 0.0_dp
      real(dp) :: cosine = 1.0_dp
   end type cyclic_encoding_type

contains

   pure elemental logical function date_valid(self) result(valid)
      !! Reports whether a date has valid Gregorian components.
      class(date_type), intent(in) :: self !! Date to validate.

      valid = self%month >= 1 .and. self%month <= 12
      if (valid) valid = self%day >= 1 .and. self%day <= month_length(self%year, self%month)
   end function date_valid

   pure elemental logical function datetime_valid(self) result(valid)
      !! Reports whether all civil clock and offset components are valid.
      class(datetime_type), intent(in) :: self !! Date-time to validate.

      valid = self%date%valid() .and. self%hour >= 0 .and. self%hour <= 23
      valid = valid .and. self%minute >= 0 .and. self%minute <= 59
      valid = valid .and. self%second >= 0.0_dp .and. self%second < 60.0_dp
      valid = valid .and. abs(self%offset_minutes) <= 24 * 60
   end function datetime_valid

   pure elemental function date_iso8601(self) result(text)
      !! Formats a valid date as YYYY-MM-DD.
      class(date_type), intent(in) :: self !! Date to format.
      character(len=10) :: text

      if (.not. self%valid()) then
         text = "NA        "
      else
         write (text, '(i4.4,"-",i2.2,"-",i2.2)') self%year, self%month, self%day
      end if
   end function date_iso8601

   pure elemental function datetime_iso8601(self) result(text)
      !! Formats a fixed-offset date-time in ISO-8601 form to microsecond precision.
      class(datetime_type), intent(in) :: self !! Date-time to format.
      character(len=32) :: text
      character(len=1) :: sign
      integer :: offset_hour, offset_minute, whole_second, microsecond

      if (.not. self%valid()) then
         text = "NA"
         return
      end if
      whole_second = int(self%second)
      microsecond = min(999999, nint((self%second - real(whole_second, dp)) * 1.0e6_dp))
      if (self%offset_minutes == 0) then
         write (text, '(i4.4,"-",i2.2,"-",i2.2,"T",i2.2,":",i2.2,":",i2.2,".",i6.6,"Z")') &
            self%date%year, self%date%month, self%date%day, self%hour, self%minute, whole_second, microsecond
      else
         sign = merge("+", "-", self%offset_minutes >= 0)
         offset_hour = abs(self%offset_minutes) / 60
         offset_minute = mod(abs(self%offset_minutes), 60)
         write (text, '(i4.4,"-",i2.2,"-",i2.2,"T",i2.2,":",i2.2,":",i2.2,".",i6.6,a,i2.2,":",i2.2)') &
            self%date%year, self%date%month, self%date%day, self%hour, self%minute, whole_second, microsecond, &
            sign, offset_hour, offset_minute
      end if
   end function datetime_iso8601

   pure elemental integer function month_length(year, month) result(number_of_days)
      !! Returns a Gregorian month length or zero for an invalid month.
      integer, intent(in) :: year  !! Gregorian year.
      integer, intent(in) :: month !! Month number.
      integer, parameter :: lengths(12) = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

      if (month < 1 .or. month > 12) then
         number_of_days = 0
      else
         number_of_days = lengths(month)
         if (month == 2 .and. leap(year)) number_of_days = 29
      end if
   end function month_length

   pure elemental logical function leap(year) result(is_leap)
      !! Implements the proleptic-Gregorian leap-year rule.
      integer, intent(in) :: year !! Gregorian year.

      is_leap = mod(year, 4) == 0 .and. (mod(year, 100) /= 0 .or. mod(year, 400) == 0)
   end function leap

end module lubridate_types
