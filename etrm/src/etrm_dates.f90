! SPDX-License-Identifier: MIT
! Derived from etrm 1.0.2, Copyright (c) 2021 etrm authors.
module etrm_dates
   use iso_fortran_env, only : int64
   implicit none
   private

   public :: civil_to_day, day_offset

contains

   pure integer(int64) function civil_to_day(year, month, day) result(serial)
      integer, intent(in) :: year, month, day
      integer(int64) :: y, m, d, era, yoe, mp, doy, doe

      y = int(year, int64)
      m = int(month, int64)
      d = int(day, int64)
      if (m <= 2_int64) y = y - 1_int64
      era = y / 400_int64
      yoe = y - era * 400_int64
      if (m > 2_int64) then
         mp = m - 3_int64
      else
         mp = m + 9_int64
      end if
      doy = (153_int64 * mp + 2_int64) / 5_int64 + d - 1_int64
      doe = yoe * 365_int64 + yoe / 4_int64 - yoe / 100_int64 + doy
      serial = era * 146097_int64 + doe - 719468_int64
   end function civil_to_day

   pure integer function day_offset(year, month, day, base_year, base_month, base_day) result(offset)
      integer, intent(in) :: year, month, day
      integer, intent(in) :: base_year, base_month, base_day
      offset = int(civil_to_day(year, month, day) - civil_to_day(base_year, base_month, base_day))
   end function day_offset

end module etrm_dates
