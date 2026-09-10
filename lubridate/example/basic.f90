! SPDX-License-Identifier: MIT
! SPDX-FileComment: Basic example for the Fortran lubridate translation.
program basic
   !! Demonstrates parsing, calendar arithmetic, and fixed-offset conversion.
   use lubridate
   implicit none

   type(date_type) :: billing_date
   type(datetime_type) :: meeting, remote_view

   billing_date = ymd('2024-01-31') + months(1)
   meeting = ymd_hms('2024-06-19 14:30:00')
   remote_view = with_tz(meeting, -300)

   print '(a)', 'Rolled billing date: ' // billing_date%iso8601()
   print '(a)', 'UTC meeting:         ' // trim(meeting%iso8601())
   print '(a)', 'UTC-05 view:         ' // trim(remote_view%iso8601())
end program basic
