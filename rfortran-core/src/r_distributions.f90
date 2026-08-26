! SPDX-License-Identifier: MIT
module r_distributions
   use, intrinsic :: ieee_arithmetic, only : ieee_negative_inf, ieee_positive_inf
   use, intrinsic :: ieee_arithmetic, only : ieee_quiet_nan, ieee_value
   use r_kinds, only : dp, r_pi
   use r_transforms, only : r_expm1, r_log1p
   implicit none
   private

   real(dp), parameter :: sqrt_two = sqrt(2.0_dp)
   real(dp), parameter :: log_sqrt_two_pi = 0.5_dp*log(2.0_dp*r_pi)

   public :: r_dnorm, r_pnorm, r_qnorm

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

   pure elemental real(dp) function r_qnorm(probability, location, scale, lower_tail, log_probability) result(value)
      real(dp), intent(in) :: probability
      real(dp), intent(in), optional :: location, scale
      logical, intent(in), optional :: lower_tail, log_probability
      real(dp) :: center, log_lower, log_upper, sigma, z
      logical :: lower, probability_is_log

      center = 0.0_dp
      if (present(location)) center = location
      sigma = 1.0_dp
      if (present(scale)) sigma = scale
      lower = .true.
      if (present(lower_tail)) lower = lower_tail
      probability_is_log = .false.
      if (present(log_probability)) probability_is_log = log_probability
      if (sigma < 0.0_dp) then
         value = ieee_value(0.0_dp, ieee_quiet_nan)
         return
      end if

      call probability_logs(probability, lower, probability_is_log, log_lower, log_upper)
      if (log_lower /= log_lower .or. log_upper /= log_upper) then
         value = ieee_value(0.0_dp, ieee_quiet_nan)
         return
      end if
      if (sigma == 0.0_dp) then
         value = center
         return
      end if
      if (log_lower == ieee_value(0.0_dp, ieee_negative_inf)) then
         value = ieee_value(0.0_dp, ieee_negative_inf)
         return
      end if
      if (log_upper == ieee_value(0.0_dp, ieee_negative_inf)) then
         value = ieee_value(0.0_dp, ieee_positive_inf)
         return
      end if

      z = inverse_standard_normal(log_lower, log_upper)
      value = center + sigma*z
   end function r_qnorm

   pure elemental subroutine probability_logs(probability, lower_tail, log_probability, log_lower, log_upper)
      real(dp), intent(in) :: probability
      logical, intent(in) :: lower_tail, log_probability
      real(dp), intent(out) :: log_lower, log_upper
      real(dp) :: log_complement

      if (log_probability) then
         if (probability > 0.0_dp .or. probability /= probability) then
            log_lower = ieee_value(0.0_dp, ieee_quiet_nan)
            log_upper = log_lower
            return
         end if
         if (probability == 0.0_dp) then
            log_complement = ieee_value(0.0_dp, ieee_negative_inf)
         else if (probability == ieee_value(0.0_dp, ieee_negative_inf)) then
            log_complement = 0.0_dp
         else
            log_complement = log(-r_expm1(probability))
         end if
         if (lower_tail) then
            log_lower = probability
            log_upper = log_complement
         else
            log_lower = log_complement
            log_upper = probability
         end if
      else
         if (probability < 0.0_dp .or. probability > 1.0_dp .or. probability /= probability) then
            log_lower = ieee_value(0.0_dp, ieee_quiet_nan)
            log_upper = log_lower
            return
         end if
         if (probability == 0.0_dp) then
            if (lower_tail) then
               log_lower = ieee_value(0.0_dp, ieee_negative_inf)
               log_upper = 0.0_dp
            else
               log_lower = 0.0_dp
               log_upper = ieee_value(0.0_dp, ieee_negative_inf)
            end if
         else if (probability == 1.0_dp) then
            if (lower_tail) then
               log_lower = 0.0_dp
               log_upper = ieee_value(0.0_dp, ieee_negative_inf)
            else
               log_lower = ieee_value(0.0_dp, ieee_negative_inf)
               log_upper = 0.0_dp
            end if
         else if (lower_tail) then
            log_lower = log(probability)
            log_upper = r_log1p(-probability)
         else
            log_lower = r_log1p(-probability)
            log_upper = log(probability)
         end if
      end if
   end subroutine probability_logs

   pure elemental real(dp) function inverse_standard_normal(log_lower, log_upper) result(value)
      real(dp), intent(in) :: log_lower, log_upper
      real(dp), parameter :: a(6) = [-3.969683028665376e1_dp, 2.209460984245205e2_dp, &
         -2.759285104469687e2_dp, 1.383577518672690e2_dp, -3.066479806614716e1_dp, &
         2.506628277459239_dp]
      real(dp), parameter :: b(5) = [-5.447609879822406e1_dp, 1.615858368580409e2_dp, &
         -1.556989798598866e2_dp, 6.680131188771972e1_dp, -1.328068155288572e1_dp]
      real(dp), parameter :: c(6) = [-7.784894002430293e-3_dp, -3.223964580411365e-1_dp, &
         -2.400758277161838_dp, -2.549732539343734_dp, 4.374664141464968_dp, &
         2.938163982698783_dp]
      real(dp), parameter :: d(4) = [7.784695709041462e-3_dp, 3.224671290700398e-1_dp, &
         2.445134137142996_dp, 3.754408661907416_dp]
      real(dp), parameter :: log_plow = log(0.02425_dp)
      real(dp) :: q, r

      if (log_lower < log_plow) then
         q = sqrt(-2.0_dp*log_lower)
         value = (((((c(1)*q + c(2))*q + c(3))*q + c(4))*q + c(5))*q + c(6)) / &
            ((((d(1)*q + d(2))*q + d(3))*q + d(4))*q + 1.0_dp)
      else if (log_upper < log_plow) then
         q = sqrt(-2.0_dp*log_upper)
         value = -(((((c(1)*q + c(2))*q + c(3))*q + c(4))*q + c(5))*q + c(6)) / &
            ((((d(1)*q + d(2))*q + d(3))*q + d(4))*q + 1.0_dp)
      else
         q = exp(log_lower) - 0.5_dp
         r = q*q
         value = (((((a(1)*r + a(2))*r + a(3))*r + a(4))*r + a(5))*r + a(6))*q / &
            (((((b(1)*r + b(2))*r + b(3))*r + b(4))*r + b(5))*r + 1.0_dp)
      end if
      call refine_normal_quantile(value, log_lower, log_upper)
   end function inverse_standard_normal

   pure elemental subroutine refine_normal_quantile(value, log_lower, log_upper)
      real(dp), intent(inout) :: value
      real(dp), intent(in) :: log_lower, log_upper
      real(dp) :: error, log_estimate, ratio, u
      integer :: i

      do i = 1, 2
         if (value < -35.0_dp) then
            log_estimate = log_normal_upper_tail(-value)
            ratio = exp(-0.5_dp*value*value - log_sqrt_two_pi - log_estimate)
            value = value - (log_estimate - log_lower)/ratio
         else if (value > 35.0_dp) then
            log_estimate = log_normal_upper_tail(value)
            ratio = exp(-0.5_dp*value*value - log_sqrt_two_pi - log_estimate)
            value = value + (log_estimate - log_upper)/ratio
         else if (value <= 0.0_dp) then
            error = 0.5_dp*erfc(-value/sqrt_two) - exp(log_lower)
            u = error/exp(-0.5_dp*value*value - log_sqrt_two_pi)
            value = value - u/(1.0_dp + 0.5_dp*value*u)
         else
            error = 0.5_dp*erfc(value/sqrt_two) - exp(log_upper)
            u = error/exp(-0.5_dp*value*value - log_sqrt_two_pi)
            value = value + u/(1.0_dp - 0.5_dp*value*u)
         end if
      end do
   end subroutine refine_normal_quantile

   pure elemental real(dp) function log_normal_upper_tail(x) result(value)
      real(dp), intent(in) :: x
      real(dp) :: inverse_square, series

      inverse_square = 1.0_dp/(x*x)
      series = 1.0_dp - inverse_square*(1.0_dp - inverse_square*(3.0_dp - &
         inverse_square*(15.0_dp - inverse_square*(105.0_dp - inverse_square*945.0_dp))))
      value = -0.5_dp*x*x - log(x) - log_sqrt_two_pi + log(series)
   end function log_normal_upper_tail

end module r_distributions
