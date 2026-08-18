! SPDX-License-Identifier: GPL-2.0-or-later
module lps_special
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
   use lps_kinds, only : dp
   implicit none
   private
   public :: lambert_wm1
   public :: expm1_dp, log1p_dp
   public :: log_expm1_pos
   public :: log_one_plus_p_expm1

contains

   pure elemental function expm1_dp(x) result(y)
      real(dp), intent(in) :: x
      real(dp) :: y, ax

      ax = abs(x)
      if (ax < 1.0e-5_dp) then
         y = x * (1.0_dp + x * (0.5_dp + x * (1.0_dp/6.0_dp + &
             x * (1.0_dp/24.0_dp + x * (1.0_dp/120.0_dp + &
             x / 720.0_dp)))))
      else
         y = exp(x) - 1.0_dp
      end if
   end function expm1_dp

   pure elemental function log1p_dp(x) result(y)
      real(dp), intent(in) :: x
      real(dp) :: y

      if (abs(x) < 1.0e-4_dp) then
         y = x * (1.0_dp + x * (-0.5_dp + x * (1.0_dp/3.0_dp + &
             x * (-0.25_dp + x * (0.2_dp - x/6.0_dp)))))
      else
         y = log(1.0_dp + x)
      end if
   end function log1p_dp

   pure elemental function lambert_wm1(z) result(w)
      ! Real W_{-1}(z) branch for -1/e <= z < 0.
      real(dp), intent(in) :: z
      real(dp) :: w
      real(dp) :: q, ew, f, denom, step, one_over_e
      integer :: iter

      one_over_e = exp(-1.0_dp)

      if (z > 0.0_dp .or. z < -one_over_e) then
         w = ieee_value(0.0_dp, ieee_quiet_nan)
         return
      end if
      if (abs(z + one_over_e) <= 4.0_dp*epsilon(1.0_dp)) then
         w = -1.0_dp
         return
      end if
      if (z >= 0.0_dp) then
         w = -huge(1.0_dp)
         return
      end if

      ! Series about the branch point is accurate and well conditioned there.
      if (z < -0.25_dp) then
         q = sqrt(max(0.0_dp, 2.0_dp * (1.0_dp + exp(1.0_dp) * z)))
         w = -1.0_dp - q - q*q / 3.0_dp - 11.0_dp*q*q*q / 72.0_dp
      else
         ! Asymptotic starting value for z -> 0-.
         w = log(-z)
         w = w - log(-w)
         if (w > -1.0_dp) w = -1.1_dp
      end if

      do iter = 1, 50
         ew = exp(w)
         f = w * ew - z
         denom = ew * (w + 1.0_dp)
         if (abs(w + 1.0_dp) > sqrt(epsilon(1.0_dp))) then
            denom = denom - (w + 2.0_dp) * f / (2.0_dp * (w + 1.0_dp))
         end if
         if (abs(denom) <= tiny(1.0_dp)) exit
         step = f / denom
         if (w - step >= -1.0_dp) step = 0.5_dp * (w + 1.0_dp)
         w = w - step
         if (abs(step) <= 8.0_dp * epsilon(1.0_dp) * max(1.0_dp, abs(w))) exit
      end do
   end function lambert_wm1

   pure elemental function log_expm1_pos(x) result(y)
      ! log(exp(x)-1) for x > 0, avoiding overflow/cancellation.
      real(dp), intent(in) :: x
      real(dp) :: y

      if (x < log(2.0_dp)) then
         y = log(expm1_dp(x))
      else
         y = x + log1p_dp(-exp(-x))
      end if
   end function log_expm1_pos

   pure elemental function log_one_plus_p_expm1(p, x) result(y)
      ! log(1 + p*(exp(x)-1)), 0 <= p <= 1, x >= 0.
      real(dp), intent(in) :: p, x
      real(dp) :: y

      if (p <= 0.0_dp) then
         y = 0.0_dp
      else if (p >= 1.0_dp) then
         y = x
      else if (x < 40.0_dp) then
         y = log1p_dp(p * expm1_dp(x))
      else
         y = x + log(p + (1.0_dp - p) * exp(-x))
      end if
   end function log_one_plus_p_expm1

end module lps_special
