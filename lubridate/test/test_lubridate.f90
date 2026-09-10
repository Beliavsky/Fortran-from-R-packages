! SPDX-License-Identifier: MIT
! SPDX-FileComment: Deterministic tests for the Fortran lubridate translation.
program test_lubridate
   !! Verifies deterministic calendar, parsing, span, rounding, and zone behavior.
   use, intrinsic :: iso_fortran_env, only: int64
   use lubridate
   implicit none

   integer :: failures
   type(cyclic_encoding_type) :: encoded
   type(date_type) :: date_value
   type(datetime_type) :: datetime_value, converted
   type(duration_type) :: duration_value
   type(interval_type) :: interval_value
   type(period_type) :: period_value

   failures = 0

   call check(leap_year(2000) .and. .not. leap_year(1900), 'Gregorian leap years', failures)
   call check(days_in_month(2024, 2) == 29, 'February in a leap year', failures)
   call check(date_to_day_number(make_date(1970, 1, 1)) == 0_int64, 'Unix epoch date', failures)
   call check(date_from_day_number(-1_int64) == make_date(1969, 12, 31), 'pre-epoch date', failures)
   call check(yday(make_date(2024, 3, 1)) == 61, 'leap-year day number', failures)
   call check(wday(make_date(1970, 1, 1)) == 5, 'Sunday-based weekday', failures)
   call check(wday(make_date(1970, 1, 1), 2) == 4, 'Monday-based weekday', failures)
   call check(isoyear(make_date(2021, 1, 1)) == 2020, 'ISO boundary year', failures)
   call check(isoweek(make_date(2021, 1, 1)) == 53, 'ISO boundary week', failures)
   call check(quarter(make_date(2024, 11, 2)) == 4, 'calendar quarter', failures)

   call check(ymd('2024-02-29') == make_date(2024, 2, 29), 'ymd parser', failures)
   call check(mdy('02/29/2024') == make_date(2024, 2, 29), 'mdy parser', failures)
   call check(dmy('29.02.2024') == make_date(2024, 2, 29), 'dmy parser', failures)
   call check(ymd('20240229') == make_date(2024, 2, 29), 'compact ymd parser', failures)
   call check(yq('2024 Q3') == make_date(2024, 7, 1), 'year-quarter parser', failures)
   datetime_value = ymd_hms('2024-02-29 23:15:07')
   call check(datetime_value%date == make_date(2024, 2, 29), 'date-time parser date', failures)
   call check(datetime_value%hour == 23 .and. datetime_value%minute == 15, 'date-time parser clock', failures)
   datetime_value = ymd_hms('20240229231507')
   call check(datetime_value%hour == 23 .and. abs(datetime_value%second - 7.0_dp) < 1.0e-12_dp, &
              'compact date-time parser', failures)
   datetime_value = ymd_hms('2024-02-29 23:15:07.25')
   call check_close(datetime_value%second, 7.25_dp, 1.0e-12_dp, 'fractional-second parser', failures)

   date_value = make_date(2023, 1, 31) + months(1)
   call check(date_value == make_date(2023, 2, 28), 'month rollback', failures)
   date_value = make_date(2024, 1, 31) + months(1)
   call check(date_value == make_date(2024, 2, 29), 'leap month rollback', failures)
   call check(rollback(make_date(2024, 3, 15)) == make_date(2024, 2, 29), 'rollback', failures)
   call check(rollforward(make_date(2024, 12, 2)) == make_date(2025, 1, 1), 'rollforward', failures)

   duration_value = dhours(1.5_dp)
   call check_close(duration_value%seconds, 5400.0_dp, 1.0e-12_dp, 'duration constructor', failures)
   period_value = period(years=1, months=2, days=3, hours=4)
   call check(period_value%years == 1 .and. period_value%months == 2, 'period constructor', failures)
   period_value = seconds_to_period(90061.5_dp)
   call check(period_value%days == 1, 'seconds-to-period days', failures)

   datetime_value = make_datetime(1970, 1, 1, 0, 0, 0.0_dp, 0)
   call check_close(datetime_to_epoch(datetime_value), 0.0_dp, 1.0e-12_dp, 'Unix epoch instant', failures)
   converted = with_tz(datetime_value, -300)
   call check(converted%date == make_date(1969, 12, 31), 'fixed-offset date conversion', failures)
   call check(converted%hour == 19, 'fixed-offset hour conversion', failures)
   call check_close(datetime_to_epoch(converted), 0.0_dp, 1.0e-12_dp, 'zone conversion preserves instant', failures)
   converted = force_tz(datetime_value, 60)
   call check_close(datetime_to_epoch(converted), -3600.0_dp, 1.0e-12_dp, 'force zone preserves clock', failures)

   interval_value = interval(datetime_value, datetime_value + dhours(2.0_dp))
   call check_close(int_length(interval_value), 7200.0_dp, 1.0e-12_dp, 'interval length', failures)
   call check(within(datetime_value + dminutes(30.0_dp), interval_value), 'interval membership', failures)
   call check_close(time_length(interval_value, 'hours'), 2.0_dp, 1.0e-12_dp, 'interval units', failures)

   datetime_value = make_datetime(2024, 6, 19, 12, 34, 40.0_dp)
   converted = floor_date(datetime_value, 'hour')
   call check(converted%hour == 12 .and. converted%minute == 0, 'floor hour', failures)
   converted = round_date(datetime_value, 'minute')
   call check(converted%minute == 35 .and. abs(converted%second) < 1.0e-12_dp, 'round minute', failures)
   call check(floor_date(make_date(2024, 6, 19), 'month') == make_date(2024, 6, 1), 'floor month', failures)

   encoded = cyclic_encoding(3.0_dp, 12.0_dp)
   call check_close(encoded%sine, 1.0_dp, 1.0e-12_dp, 'cyclic sine', failures)
   call check_close(encoded%cosine, 0.0_dp, 1.0e-12_dp, 'cyclic cosine', failures)

   if (failures > 0) error stop 'lubridate tests failed'
   print '(a)', 'All lubridate tests passed.'

contains

   subroutine check(condition, label, failure_count)
      !! Records a failed logical assertion.
      logical, intent(in) :: condition      !! Assertion result.
      character(len=*), intent(in) :: label !! Assertion label.
      integer, intent(inout) :: failure_count !! Accumulated failure count.
      if (.not. condition) then
         print '(a)', 'FAIL: ' // label
         failure_count = failure_count + 1
      end if
   end subroutine check

   subroutine check_close(actual, expected, tolerance, label, failure_count)
      !! Records a failed approximate numeric assertion.
      real(dp), intent(in) :: actual         !! Computed value.
      real(dp), intent(in) :: expected       !! Reference value.
      real(dp), intent(in) :: tolerance      !! Absolute tolerance.
      character(len=*), intent(in) :: label  !! Assertion label.
      integer, intent(inout) :: failure_count !! Accumulated failure count.
      call check(abs(actual - expected) <= tolerance, label, failure_count)
   end subroutine check_close

end program test_lubridate
