! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) 2013 Bruno Remillard
! Modern Fortran translation copyright (C) 2026 OpenAI
module opthedging_interpolation
   use opthedging_kinds, only : dp
   implicit none
   private

   public :: interpolation1d
   public :: linear_interpolate_uniform

contains

   pure function linear_interpolate_uniform(x, values, min_x, max_x) result(y)
      real(dp), intent(in) :: x
      real(dp), intent(in) :: values(:)
      real(dp), intent(in) :: min_x
      real(dp), intent(in) :: max_x
      real(dp) :: y

      integer :: j
      integer :: n
      real(dp) :: dx
      real(dp) :: position
      real(dp) :: weight

      n = size(values)
      if (n == 0) then
         y = 0.0_dp
         return
      end if
      if (n == 1 .or. abs(max_x - min_x) <= tiny(1.0_dp)) then
         y = values(1)
         return
      end if

      dx = (max_x - min_x) / real(n - 1, dp)
      if (x <= min_x) then
         y = values(1) + (x - min_x) * (values(2) - values(1)) / dx
      else if (x >= max_x) then
         y = values(n) + (x - max_x) * (values(n) - values(n - 1)) / dx
      else
         position = (x - min_x) / dx
         j = min(int(floor(position)) + 1, n - 1)
         weight = position - real(j - 1, dp)
         y = (1.0_dp - weight) * values(j) + weight * values(j + 1)
      end if
   end function linear_interpolate_uniform

   pure function interpolation1d(x, values, min_s, max_s) result(y)
      real(dp), intent(in) :: x
      real(dp), intent(in) :: values(:)
      real(dp), intent(in) :: min_s
      real(dp), intent(in) :: max_s
      real(dp) :: y

      y = linear_interpolate_uniform(x, values, min_s, max_s)
   end function interpolation1d

end module opthedging_interpolation
