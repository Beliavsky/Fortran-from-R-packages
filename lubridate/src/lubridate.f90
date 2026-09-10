! SPDX-License-Identifier: MIT
! SPDX-FileComment: Public facade for the Fortran lubridate translation.
module lubridate
   !! Provides a modern Fortran subset of lubridate's date-time API.
   use lubridate_types
   use lubridate_calendar
   use lubridate_parse
   use lubridate_spans
   use lubridate_round
   implicit none
   public
end module lubridate
