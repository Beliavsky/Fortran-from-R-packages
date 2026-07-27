! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) 2013 Bruno Remillard
! Modern Fortran translation copyright (C) 2026 OpenAI
module opthedging_statistics
   use opthedging_kinds, only : dp
   implicit none
   private

   public :: mean_square
   public :: mean_value

contains

   pure function mean_value(x) result(value)
      real(dp), intent(in) :: x(:)
      real(dp) :: value

      if (size(x) == 0) then
         value = 0.0_dp
      else
         value = sum(x) / real(size(x), dp)
      end if
   end function mean_value

   pure function mean_square(x) result(value)
      real(dp), intent(in) :: x(:)
      real(dp) :: value

      if (size(x) == 0) then
         value = 0.0_dp
      else
         value = dot_product(x, x) / real(size(x), dp)
      end if
   end function mean_square

end module opthedging_statistics
