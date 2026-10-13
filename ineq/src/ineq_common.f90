! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
module ineq_common
   !! Internal input handling for the ineq translation.
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_is_finite, ieee_value, ieee_quiet_nan
   use r_kinds, only : dp
   use r_ordering, only : r_sort_values_in_place
   implicit none
   private
   public :: prepare_sample, normalize_sample, dp, nan_value
contains
   pure real(dp) function nan_value() result(value)
      value = ieee_value(0.0_dp, ieee_quiet_nan)
   end function nan_value

   pure subroutine prepare_sample(x, values, valid, na_rm, sort_values)
      real(dp), intent(in) :: x(:) !! Observations; NaNs represent R NA and NaN.
      real(dp), allocatable, intent(out) :: values(:) !! Retained observations, optionally sorted ascending.
      logical, intent(out) :: valid !! False for retained NaNs or infinite observations.
      logical, intent(in), optional :: na_rm !! Remove NaNs when true; defaults to true.
      logical, intent(in), optional :: sort_values !! Sort retained values when true; defaults to false.
      logical :: remove_na
      remove_na = .true.
      if (present(na_rm)) remove_na = na_rm
      values = pack(x, .not. ieee_is_nan(x))
      valid = all(ieee_is_finite(values))
      if (.not. remove_na .and. any(ieee_is_nan(x))) valid = .false.
      if (present(sort_values)) then
         if (sort_values .and. valid) call r_sort_values_in_place(values)
      end if
   end subroutine prepare_sample

   pure subroutine normalize_sample(values, scale)
      real(dp), intent(inout) :: values(:) !! Nonempty finite observations, divided by their largest absolute value.
      real(dp), intent(out) :: scale !! Positive divisor used; one for an all-zero sample.
      scale = maxval(abs(values))
      if (scale == 0.0_dp) scale = 1.0_dp
      values = values / scale
   end subroutine normalize_sample
end module ineq_common
