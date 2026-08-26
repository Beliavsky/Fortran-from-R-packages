! SPDX-License-Identifier: MIT
module r_robust
   use, intrinsic :: ieee_arithmetic, only : ieee_quiet_nan, ieee_value
   use r_kinds, only : dp
   use r_quantiles, only : r_median
   implicit none
   private

   public :: r_mad

contains

   pure real(dp) function r_mad(x, center, constant, na_rm) result(value)
      real(dp), intent(in) :: x(:)
      real(dp), intent(in), optional :: center, constant
      logical, intent(in), optional :: na_rm
      real(dp), allocatable :: deviations(:)
      real(dp) :: center_value, multiplier

      multiplier = 1.4826_dp
      if (present(constant)) multiplier = constant
      if (present(center)) then
         center_value = center
      else
         center_value = r_median(x, na_rm)
      end if
      if (center_value /= center_value .or. multiplier /= multiplier) then
         value = ieee_value(0.0_dp, ieee_quiet_nan)
         return
      end if
      allocate(deviations(size(x)))
      deviations = abs(x - center_value)
      value = multiplier*r_median(deviations, na_rm)
   end function r_mad

end module r_robust
