! SPDX-License-Identifier: MIT
module r_transforms
   use, intrinsic :: ieee_arithmetic, only : ieee_negative_inf, ieee_positive_inf, &
      ieee_quiet_nan, ieee_value
   use r_kinds, only : dp
   implicit none
   private

   public :: r_expm1, r_log1mexp, r_log1p, r_log1pexp, r_logistic, r_logit

contains

   pure elemental real(dp) function r_log1p(x) result(value)
      real(dp), intent(in) :: x
      real(dp) :: term
      integer :: k

      if (x < -1.0_dp) then
         value = ieee_value(0.0_dp, ieee_quiet_nan)
      else if (x == -1.0_dp) then
         value = ieee_value(0.0_dp, ieee_negative_inf)
      else if (abs(x) < 1.0e-4_dp) then
         value = x
         term = x
         do k = 2, 50
            term = -term*x*real(k - 1, dp)/real(k, dp)
            value = value + term
            if (abs(term) <= epsilon(1.0_dp)*max(abs(value), tiny(1.0_dp))) exit
         end do
      else
         value = log(1.0_dp + x)
      end if
   end function r_log1p

   pure elemental real(dp) function r_expm1(x) result(value)
      real(dp), intent(in) :: x
      real(dp) :: term
      integer :: k

      if (abs(x) < 1.0e-5_dp) then
         value = x
         term = x
         do k = 2, 50
            term = term*x/real(k, dp)
            value = value + term
            if (abs(term) <= epsilon(1.0_dp)*max(abs(value), tiny(1.0_dp))) exit
         end do
      else
         value = exp(x) - 1.0_dp
      end if
   end function r_expm1

   pure elemental real(dp) function r_log1mexp(x) result(value)
      real(dp), intent(in) :: x
      real(dp), parameter :: log_half = -0.693147180559945309417232121458176568_dp

      if (x > 0.0_dp) then
         value = ieee_value(0.0_dp, ieee_quiet_nan)
      else if (x == 0.0_dp) then
         value = ieee_value(0.0_dp, ieee_negative_inf)
      else if (x < log_half) then
         value = r_log1p(-exp(x))
      else
         value = log(-r_expm1(x))
      end if
   end function r_log1mexp

   pure elemental real(dp) function r_log1pexp(x) result(value)
      real(dp), intent(in) :: x

      if (x > 0.0_dp) then
         value = x + r_log1p(exp(-x))
      else
         value = r_log1p(exp(x))
      end if
   end function r_log1pexp

   pure elemental real(dp) function r_logistic(x) result(value)
      real(dp), intent(in) :: x
      real(dp) :: exponential

      if (x >= 0.0_dp) then
         value = 1.0_dp/(1.0_dp + exp(-x))
      else
         exponential = exp(x)
         value = exponential/(1.0_dp + exponential)
      end if
   end function r_logistic

   pure elemental real(dp) function r_logit(p) result(value)
      real(dp), intent(in) :: p

      if (p < 0.0_dp .or. p > 1.0_dp) then
         value = ieee_value(0.0_dp, ieee_quiet_nan)
      else if (p == 0.0_dp) then
         value = ieee_value(0.0_dp, ieee_negative_inf)
      else if (p == 1.0_dp) then
         value = ieee_value(0.0_dp, ieee_positive_inf)
      else
         value = log(p) - r_log1p(-p)
      end if
   end function r_logit

end module r_transforms
