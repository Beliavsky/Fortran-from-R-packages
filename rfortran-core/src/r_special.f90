! SPDX-License-Identifier: MIT
module r_special
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_quiet_nan, ieee_value
   use r_kinds, only : dp
   implicit none
   private

   public :: r_digamma, r_trigamma, r_log_beta

contains

   pure elemental real(dp) function r_digamma(x) result(value)
      real(dp), intent(in) :: x
      real(dp) :: shifted, inverse, inverse_squared

      if (x <= 0.0_dp .or. .not. ieee_is_finite(x)) then
         value = ieee_value(0.0_dp, ieee_quiet_nan)
         return
      end if

      value = 0.0_dp
      shifted = x
      do while (shifted < 8.0_dp)
         value = value - 1.0_dp/shifted
         shifted = shifted + 1.0_dp
      end do
      inverse = 1.0_dp/shifted
      inverse_squared = inverse*inverse
      value = value + log(shifted) - 0.5_dp*inverse - inverse_squared*(1.0_dp/12.0_dp - &
         inverse_squared*(1.0_dp/120.0_dp - inverse_squared*(1.0_dp/252.0_dp - &
         inverse_squared*(1.0_dp/240.0_dp - inverse_squared*(5.0_dp/660.0_dp)))))
   end function r_digamma

   pure elemental real(dp) function r_trigamma(x) result(value)
      real(dp), intent(in) :: x
      real(dp) :: shifted, inverse, inverse_squared

      if (x <= 0.0_dp .or. .not. ieee_is_finite(x)) then
         value = ieee_value(0.0_dp, ieee_quiet_nan)
         return
      end if

      value = 0.0_dp
      shifted = x
      do while (shifted < 8.0_dp)
         value = value + 1.0_dp/(shifted*shifted)
         shifted = shifted + 1.0_dp
      end do
      inverse = 1.0_dp/shifted
      inverse_squared = inverse*inverse
      value = value + inverse + 0.5_dp*inverse_squared + inverse*inverse_squared/6.0_dp - &
         inverse*inverse_squared**2/30.0_dp + inverse*inverse_squared**3/42.0_dp - &
         inverse*inverse_squared**4/30.0_dp + 5.0_dp*inverse*inverse_squared**5/66.0_dp
   end function r_trigamma

   pure elemental real(dp) function r_log_beta(a, b) result(value)
      real(dp), intent(in) :: a, b

      if (a <= 0.0_dp .or. b <= 0.0_dp .or. .not. ieee_is_finite(a) .or. &
          .not. ieee_is_finite(b)) then
         value = ieee_value(0.0_dp, ieee_quiet_nan)
      else
         value = log_gamma(a) + log_gamma(b) - log_gamma(a + b)
      end if
   end function r_log_beta

end module r_special
