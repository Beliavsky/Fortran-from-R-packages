! SPDX-License-Identifier: GPL-3.0-only
module scs_linalg
   use scs_kinds, only : dp, i4
   implicit none
   private
   public :: norm_2, norm_inf, norm_sq, norm_diff, mean_value, safe_div_pos
contains
   pure real(dp) function norm_sq(x) result(v)
      real(dp), intent(in) :: x(:)
      v = dot_product(x, x)
   end function norm_sq

   pure real(dp) function norm_2(x) result(v)
      real(dp), intent(in) :: x(:)
      v = sqrt(max(0.0_dp, dot_product(x, x)))
   end function norm_2

   pure real(dp) function norm_inf(x) result(v)
      real(dp), intent(in) :: x(:)
      if (size(x) == 0) then
         v = 0.0_dp
      else
         v = maxval(abs(x))
      end if
   end function norm_inf

   pure real(dp) function norm_diff(x, y) result(v)
      real(dp), intent(in) :: x(:), y(:)
      v = norm_2(x - y)
   end function norm_diff

   pure real(dp) function mean_value(x) result(v)
      real(dp), intent(in) :: x(:)
      if (size(x) == 0) then
         v = 0.0_dp
      else
         v = sum(x) / real(size(x), dp)
      end if
   end function mean_value

   pure real(dp) function safe_div_pos(a, b) result(v)
      real(dp), intent(in) :: a, b
      if (b > 0.0_dp) then
         v = a / b
      else
         v = huge(1.0_dp)
      end if
   end function safe_div_pos
end module scs_linalg
