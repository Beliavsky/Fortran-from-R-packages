! SPDX-License-Identifier: MIT
module r_descriptive
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan, ieee_negative_inf, &
      ieee_positive_inf, ieee_quiet_nan, ieee_value
   use r_kinds, only : dp
   implicit none
   private

   public :: r_count_nonmissing, r_mean, r_variance, r_sd

contains

   pure integer function r_count_nonmissing(x) result(n)
      real(dp), intent(in) :: x(:)
      integer :: i

      n = 0
      do i = 1, size(x)
         if (.not. ieee_is_nan(x(i))) n = n + 1
      end do
   end function r_count_nonmissing

   pure real(dp) function r_mean(x, na_rm) result(value)
      real(dp), intent(in) :: x(:)
      logical, intent(in), optional :: na_rm
      real(dp) :: total, correction, y, updated
      integer :: i, n
      logical :: remove_na, has_positive_inf, has_negative_inf

      remove_na = .false.
      if (present(na_rm)) remove_na = na_rm
      total = 0.0_dp
      correction = 0.0_dp
      n = 0
      has_positive_inf = .false.
      has_negative_inf = .false.

      do i = 1, size(x)
         if (ieee_is_nan(x(i))) then
            if (remove_na) cycle
            value = ieee_value(0.0_dp, ieee_quiet_nan)
            return
         end if
         n = n + 1
         if (.not. ieee_is_finite(x(i))) then
            if (x(i) > 0.0_dp) has_positive_inf = .true.
            if (x(i) < 0.0_dp) has_negative_inf = .true.
            cycle
         end if
         y = x(i) - correction
         updated = total + y
         correction = (updated - total) - y
         total = updated
      end do

      if (n == 0 .or. (has_positive_inf .and. has_negative_inf)) then
         value = ieee_value(0.0_dp, ieee_quiet_nan)
      else if (has_positive_inf) then
         value = ieee_value(0.0_dp, ieee_positive_inf)
      else if (has_negative_inf) then
         value = ieee_value(0.0_dp, ieee_negative_inf)
      else
         value = total / real(n, dp)
      end if
   end function r_mean

   pure real(dp) function r_variance(x, na_rm, center, ddof) result(value)
      real(dp), intent(in) :: x(:)
      logical, intent(in), optional :: na_rm
      real(dp), intent(in), optional :: center
      integer, intent(in), optional :: ddof
      real(dp) :: mean_value, total, correction, term, updated
      integer :: degrees, i, n
      logical :: remove_na

      remove_na = .false.
      if (present(na_rm)) remove_na = na_rm
      degrees = 1
      if (present(ddof)) degrees = ddof
      if (degrees < 0) then
         value = ieee_value(0.0_dp, ieee_quiet_nan)
         return
      end if

      if (present(center)) then
         mean_value = center
      else
         mean_value = r_mean(x, remove_na)
      end if
      if (.not. ieee_is_finite(mean_value)) then
         value = ieee_value(0.0_dp, ieee_quiet_nan)
         return
      end if

      total = 0.0_dp
      correction = 0.0_dp
      n = 0
      do i = 1, size(x)
         if (ieee_is_nan(x(i))) then
            if (remove_na) cycle
            value = ieee_value(0.0_dp, ieee_quiet_nan)
            return
         end if
         if (.not. ieee_is_finite(x(i))) then
            value = ieee_value(0.0_dp, ieee_quiet_nan)
            return
         end if
         n = n + 1
         term = (x(i) - mean_value)**2 - correction
         updated = total + term
         correction = (updated - total) - term
         total = updated
      end do

      if (n <= degrees) then
         value = ieee_value(0.0_dp, ieee_quiet_nan)
      else
         value = total / real(n - degrees, dp)
      end if
   end function r_variance

   pure real(dp) function r_sd(x, na_rm, center, ddof) result(value)
      real(dp), intent(in) :: x(:)
      logical, intent(in), optional :: na_rm
      real(dp), intent(in), optional :: center
      integer, intent(in), optional :: ddof

      value = r_variance(x, na_rm, center, ddof)
      if (value >= 0.0_dp) value = sqrt(value)
   end function r_sd

end module r_descriptive
