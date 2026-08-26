! SPDX-License-Identifier: MIT
module r_distributions
   use, intrinsic :: ieee_arithmetic, only : ieee_quiet_nan, ieee_value
   use r_kinds, only : dp, r_pi
   implicit none
   private

   real(dp), parameter :: sqrt_two = sqrt(2.0_dp)
   real(dp), parameter :: log_sqrt_two_pi = 0.5_dp*log(2.0_dp*r_pi)

   public :: r_dnorm, r_pnorm

contains

   pure elemental real(dp) function r_dnorm(x, location, scale, log_density) result(value)
      real(dp), intent(in) :: x
      real(dp), intent(in), optional :: location, scale
      logical, intent(in), optional :: log_density
      real(dp) :: center, sigma, z
      logical :: return_log

      center = 0.0_dp
      if (present(location)) center = location
      sigma = 1.0_dp
      if (present(scale)) sigma = scale
      return_log = .false.
      if (present(log_density)) return_log = log_density
      if (sigma <= 0.0_dp) then
         value = ieee_value(0.0_dp, ieee_quiet_nan)
         return
      end if

      z = (x - center)/sigma
      value = -0.5_dp*z*z - log(sigma) - log_sqrt_two_pi
      if (.not. return_log) value = exp(value)
   end function r_dnorm

   pure elemental real(dp) function r_pnorm(x, location, scale, lower_tail, log_probability) result(value)
      real(dp), intent(in) :: x
      real(dp), intent(in), optional :: location, scale
      logical, intent(in), optional :: lower_tail, log_probability
      real(dp) :: center, probability, sigma, z
      logical :: lower, return_log

      center = 0.0_dp
      if (present(location)) center = location
      sigma = 1.0_dp
      if (present(scale)) sigma = scale
      lower = .true.
      if (present(lower_tail)) lower = lower_tail
      return_log = .false.
      if (present(log_probability)) return_log = log_probability
      if (sigma <= 0.0_dp) then
         value = ieee_value(0.0_dp, ieee_quiet_nan)
         return
      end if

      z = (x - center)/sigma
      if (lower) then
         probability = 0.5_dp*erfc(-z/sqrt_two)
         if (return_log .and. z < -35.0_dp) then
            value = log_normal_upper_tail(-z)
         else if (return_log) then
            value = log(probability)
         else
            value = probability
         end if
      else
         probability = 0.5_dp*erfc(z/sqrt_two)
         if (return_log .and. z > 35.0_dp) then
            value = log_normal_upper_tail(z)
         else if (return_log) then
            value = log(probability)
         else
            value = probability
         end if
      end if
   end function r_pnorm

   pure elemental real(dp) function log_normal_upper_tail(x) result(value)
      real(dp), intent(in) :: x
      real(dp) :: inverse_square, series

      inverse_square = 1.0_dp/(x*x)
      series = 1.0_dp - inverse_square*(1.0_dp - inverse_square*(3.0_dp - &
         inverse_square*(15.0_dp - inverse_square*(105.0_dp - inverse_square*945.0_dp))))
      value = -0.5_dp*x*x - log(x) - log_sqrt_two_pi + log(series)
   end function log_normal_upper_tail

end module r_distributions
