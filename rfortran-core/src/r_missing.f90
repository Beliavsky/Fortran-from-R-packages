! SPDX-License-Identifier: MIT
module r_missing
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
   use r_kinds, only : dp
   implicit none
   private

   public :: r_is_na, r_is_finite

contains

   pure elemental logical function r_is_na(x) result(is_na)
      real(dp), intent(in) :: x

      is_na = ieee_is_nan(x)
   end function r_is_na

   pure elemental logical function r_is_finite(x) result(is_finite)
      real(dp), intent(in) :: x

      is_finite = ieee_is_finite(x)
   end function r_is_finite

end module r_missing
