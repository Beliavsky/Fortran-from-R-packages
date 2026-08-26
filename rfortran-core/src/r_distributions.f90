! SPDX-License-Identifier: MIT
module r_distributions
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_negative_inf, ieee_positive_inf
   use, intrinsic :: ieee_arithmetic, only : ieee_quiet_nan, ieee_value
   use r_kinds, only : dp, r_pi
   use r_special, only : r_log_beta, r_regularized_beta
   use r_special, only : r_regularized_gamma_p, r_regularized_gamma_q
   use r_transforms, only : r_expm1, r_log1p
   implicit none
   private

   real(dp), parameter :: sqrt_two = sqrt(2.0_dp)
   real(dp), parameter :: log_sqrt_two_pi = 0.5_dp*log(2.0_dp*r_pi)

   public :: r_df, r_dchisq, r_dnorm, r_dt
   public :: r_pf, r_pchisq, r_pnorm, r_pt
   public :: r_qf, r_qchisq, r_qnorm, r_qt

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

   pure elemental real(dp) function r_dt(x, degrees_freedom, log_density) result(value)
      real(dp), intent(in) :: x, degrees_freedom
      logical, intent(in), optional :: log_density
      logical :: return_log

      return_log = .false.
      if (present(log_density)) return_log = log_density
      if (degrees_freedom /= degrees_freedom .or. degrees_freedom <= 0.0_dp) then
         value = ieee_value(0.0_dp, ieee_quiet_nan)
         return
      end if
      if (.not. ieee_is_finite(degrees_freedom)) then
         value = r_dnorm(x, log_density=return_log)
         return
      end if
      value = log_gamma(0.5_dp*(degrees_freedom + 1.0_dp)) - &
         log_gamma(0.5_dp*degrees_freedom) - 0.5_dp*log(degrees_freedom*r_pi) - &
         0.5_dp*(degrees_freedom + 1.0_dp)*r_log1p(x*x/degrees_freedom)
      if (.not. return_log) value = exp(value)
   end function r_dt

   pure elemental real(dp) function r_pt(x, degrees_freedom, lower_tail, log_probability) result(value)
      real(dp), intent(in) :: x, degrees_freedom
      logical, intent(in), optional :: lower_tail, log_probability
      real(dp) :: lower_value, small_tail, upper_value, transformed
      logical :: lower, return_log

      lower = .true.
      if (present(lower_tail)) lower = lower_tail
      return_log = .false.
      if (present(log_probability)) return_log = log_probability
      if (degrees_freedom /= degrees_freedom .or. degrees_freedom <= 0.0_dp) then
         value = ieee_value(0.0_dp, ieee_quiet_nan)
         return
      end if
      if (.not. ieee_is_finite(degrees_freedom)) then
         value = r_pnorm(x, lower_tail=lower, log_probability=return_log)
         return
      end if
      if (x /= x) then
         value = ieee_value(0.0_dp, ieee_quiet_nan)
         return
      end if
      if (x == 0.0_dp) then
         lower_value = 0.5_dp
         upper_value = 0.5_dp
      else
         transformed = degrees_freedom/(degrees_freedom + x*x)
         small_tail = 0.5_dp*r_regularized_beta(transformed,0.5_dp*degrees_freedom,0.5_dp)
         if (x < 0.0_dp) then
            lower_value = small_tail
            upper_value = 1.0_dp - small_tail
         else
            lower_value = 1.0_dp - small_tail
            upper_value = small_tail
         end if
      end if
      value = probability_result(lower_value,upper_value,lower,return_log)
   end function r_pt

   pure elemental real(dp) function r_qt(probability, degrees_freedom, lower_tail, log_probability) result(value)
      real(dp), intent(in) :: probability, degrees_freedom
      logical, intent(in), optional :: lower_tail, log_probability
      real(dp) :: high, log_lower, log_small, log_upper, low, midpoint
      logical :: lower, probability_is_log
      integer :: iteration

      lower = .true.
      if (present(lower_tail)) lower = lower_tail
      probability_is_log = .false.
      if (present(log_probability)) probability_is_log = log_probability
      if (degrees_freedom /= degrees_freedom .or. degrees_freedom <= 0.0_dp) then
         value = ieee_value(0.0_dp, ieee_quiet_nan)
         return
      end if
      if (.not. ieee_is_finite(degrees_freedom)) then
         value = r_qnorm(probability,lower_tail=lower,log_probability=probability_is_log)
         return
      end if
      call probability_logs(probability,lower,probability_is_log,log_lower,log_upper)
      if (log_lower /= log_lower .or. log_upper /= log_upper) then
         value = ieee_value(0.0_dp, ieee_quiet_nan)
         return
      end if
      if (log_lower == ieee_value(0.0_dp,ieee_negative_inf)) then
         value = ieee_value(0.0_dp,ieee_negative_inf)
         return
      else if (log_upper == ieee_value(0.0_dp,ieee_negative_inf)) then
         value = ieee_value(0.0_dp,ieee_positive_inf)
         return
      else if (log_lower == log_upper) then
         value = 0.0_dp
         return
      end if

      log_small = min(log_lower,log_upper)
      high = max(1.0_dp,abs(r_qnorm(log_small,log_probability=.true.)))
      do while (r_pt(high,degrees_freedom,lower_tail=.false.,log_probability=.true.) > log_small)
         if (high >= 0.25_dp*huge(1.0_dp)) exit
         high = 2.0_dp*high
      end do
      low = 0.0_dp
      do iteration = 1, 240
         midpoint = 0.5_dp*(low + high)
         if (r_pt(midpoint,degrees_freedom,lower_tail=.false.,log_probability=.true.) > log_small) then
            low = midpoint
         else
            high = midpoint
         end if
         if (high - low <= 8.0_dp*epsilon(1.0_dp)*max(1.0_dp,midpoint)) exit
      end do
      value = 0.5_dp*(low + high)
      if (log_lower < log_upper) value = -value
   end function r_qt

   pure elemental real(dp) function r_dchisq(x, degrees_freedom, log_density) result(value)
      real(dp), intent(in) :: x, degrees_freedom
      logical, intent(in), optional :: log_density
      real(dp) :: shape
      logical :: return_log

      return_log = .false.
      if (present(log_density)) return_log = log_density
      if (degrees_freedom <= 0.0_dp .or. .not. ieee_is_finite(degrees_freedom) .or. x /= x) then
         value = ieee_value(0.0_dp,ieee_quiet_nan)
      else if (x < 0.0_dp .or. .not. ieee_is_finite(x)) then
         value = ieee_value(0.0_dp,ieee_negative_inf)
         if (.not. return_log) value = 0.0_dp
      else if (x == 0.0_dp) then
         if (degrees_freedom < 2.0_dp) then
            value = ieee_value(0.0_dp,ieee_positive_inf)
         else if (degrees_freedom == 2.0_dp) then
            value = -log(2.0_dp)
            if (.not. return_log) value = 0.5_dp
         else
            value = ieee_value(0.0_dp,ieee_negative_inf)
            if (.not. return_log) value = 0.0_dp
         end if
      else
         shape = 0.5_dp*degrees_freedom
         value = (shape - 1.0_dp)*log(x) - 0.5_dp*x - shape*log(2.0_dp) - log_gamma(shape)
         if (.not. return_log) value = exp(value)
      end if
   end function r_dchisq

   pure elemental real(dp) function r_pchisq(x, degrees_freedom, lower_tail, log_probability) result(value)
      real(dp), intent(in) :: x, degrees_freedom
      logical, intent(in), optional :: lower_tail, log_probability
      real(dp) :: lower_value, upper_value
      logical :: lower, return_log

      lower = .true.
      if (present(lower_tail)) lower = lower_tail
      return_log = .false.
      if (present(log_probability)) return_log = log_probability
      if (degrees_freedom <= 0.0_dp .or. .not. ieee_is_finite(degrees_freedom) .or. x /= x) then
         value = ieee_value(0.0_dp,ieee_quiet_nan)
         return
      end if
      if (x <= 0.0_dp) then
         lower_value = 0.0_dp
         upper_value = 1.0_dp
      else
         lower_value = r_regularized_gamma_p(0.5_dp*degrees_freedom,0.5_dp*x)
         upper_value = r_regularized_gamma_q(0.5_dp*degrees_freedom,0.5_dp*x)
      end if
      value = probability_result(lower_value,upper_value,lower,return_log)
   end function r_pchisq

   pure elemental real(dp) function r_qchisq(probability, degrees_freedom, lower_tail, log_probability) result(value)
      real(dp), intent(in) :: probability, degrees_freedom
      logical, intent(in), optional :: lower_tail, log_probability
      real(dp) :: high, log_lower, log_upper, low, midpoint
      logical :: invert_lower, lower, probability_is_log
      integer :: iteration

      lower = .true.
      if (present(lower_tail)) lower = lower_tail
      probability_is_log = .false.
      if (present(log_probability)) probability_is_log = log_probability
      if (degrees_freedom <= 0.0_dp .or. .not. ieee_is_finite(degrees_freedom)) then
         value = ieee_value(0.0_dp,ieee_quiet_nan)
         return
      end if
      call probability_logs(probability,lower,probability_is_log,log_lower,log_upper)
      if (log_lower /= log_lower .or. log_upper /= log_upper) then
         value = ieee_value(0.0_dp,ieee_quiet_nan)
         return
      end if
      if (log_lower == ieee_value(0.0_dp,ieee_negative_inf)) then
         value = 0.0_dp
         return
      else if (log_upper == ieee_value(0.0_dp,ieee_negative_inf)) then
         value = ieee_value(0.0_dp,ieee_positive_inf)
         return
      end if

      invert_lower = log_lower <= log_upper
      low = 0.0_dp
      high = max(1.0_dp,degrees_freedom)
      if (invert_lower) then
         do while (r_pchisq(high,degrees_freedom,log_probability=.true.) < log_lower)
            if (high >= 0.25_dp*huge(1.0_dp)) exit
            high = 2.0_dp*high
         end do
      else
         do while (r_pchisq(high,degrees_freedom,lower_tail=.false.,log_probability=.true.) > log_upper)
            if (high >= 0.25_dp*huge(1.0_dp)) exit
            high = 2.0_dp*high
         end do
      end if
      do iteration = 1, 240
         midpoint = 0.5_dp*(low + high)
         if (invert_lower) then
            if (r_pchisq(midpoint,degrees_freedom,log_probability=.true.) < log_lower) then
               low = midpoint
            else
               high = midpoint
            end if
         else
            if (r_pchisq(midpoint,degrees_freedom,lower_tail=.false.,log_probability=.true.) > log_upper) then
               low = midpoint
            else
               high = midpoint
            end if
         end if
         if (high - low <= 8.0_dp*epsilon(1.0_dp)*max(1.0_dp,midpoint)) exit
      end do
      value = 0.5_dp*(low + high)
   end function r_qchisq

   pure elemental real(dp) function r_df(x, degrees_freedom1, degrees_freedom2, log_density) result(value)
      real(dp), intent(in) :: x, degrees_freedom1, degrees_freedom2
      logical, intent(in), optional :: log_density
      real(dp) :: shape1, shape2
      logical :: return_log

      return_log = .false.
      if (present(log_density)) return_log = log_density
      if (degrees_freedom1 <= 0.0_dp .or. degrees_freedom2 <= 0.0_dp .or. &
          .not. ieee_is_finite(degrees_freedom1) .or. .not. ieee_is_finite(degrees_freedom2) .or. x /= x) then
         value = ieee_value(0.0_dp,ieee_quiet_nan)
         return
      end if
      if (x < 0.0_dp .or. .not. ieee_is_finite(x)) then
         value = ieee_value(0.0_dp,ieee_negative_inf)
         if (.not. return_log) value = 0.0_dp
         return
      end if
      shape1 = 0.5_dp*degrees_freedom1
      shape2 = 0.5_dp*degrees_freedom2
      if (x == 0.0_dp) then
         if (degrees_freedom1 < 2.0_dp) then
            value = ieee_value(0.0_dp,ieee_positive_inf)
         else if (degrees_freedom1 == 2.0_dp) then
            value = shape1*log(degrees_freedom1/degrees_freedom2) - r_log_beta(shape1,shape2)
            if (.not. return_log) value = exp(value)
         else
            value = ieee_value(0.0_dp,ieee_negative_inf)
            if (.not. return_log) value = 0.0_dp
         end if
         return
      end if
      value = shape1*log(degrees_freedom1/degrees_freedom2) + (shape1 - 1.0_dp)*log(x) - &
         (shape1 + shape2)*r_log1p(degrees_freedom1*x/degrees_freedom2) - r_log_beta(shape1,shape2)
      if (.not. return_log) value = exp(value)
   end function r_df

   pure elemental real(dp) function r_pf(x, degrees_freedom1, degrees_freedom2, lower_tail, &
      log_probability) result(value)
      real(dp), intent(in) :: x, degrees_freedom1, degrees_freedom2
      logical, intent(in), optional :: lower_tail, log_probability
      real(dp) :: lower_value, transformed, upper_value
      logical :: lower, return_log

      lower = .true.
      if (present(lower_tail)) lower = lower_tail
      return_log = .false.
      if (present(log_probability)) return_log = log_probability
      if (degrees_freedom1 <= 0.0_dp .or. degrees_freedom2 <= 0.0_dp .or. &
          .not. ieee_is_finite(degrees_freedom1) .or. .not. ieee_is_finite(degrees_freedom2) .or. x /= x) then
         value = ieee_value(0.0_dp,ieee_quiet_nan)
         return
      end if
      if (x <= 0.0_dp) then
         lower_value = 0.0_dp
         upper_value = 1.0_dp
      else if (.not. ieee_is_finite(x)) then
         lower_value = 1.0_dp
         upper_value = 0.0_dp
      else
         transformed = 1.0_dp/(1.0_dp + degrees_freedom2/(degrees_freedom1*x))
         lower_value = r_regularized_beta(transformed,0.5_dp*degrees_freedom1,0.5_dp*degrees_freedom2)
         transformed = degrees_freedom2/(degrees_freedom2 + degrees_freedom1*x)
         upper_value = r_regularized_beta(transformed,0.5_dp*degrees_freedom2,0.5_dp*degrees_freedom1)
      end if
      value = probability_result(lower_value,upper_value,lower,return_log)
   end function r_pf

   pure elemental real(dp) function r_qf(probability, degrees_freedom1, degrees_freedom2, &
      lower_tail, log_probability) result(value)
      real(dp), intent(in) :: probability, degrees_freedom1, degrees_freedom2
      logical, intent(in), optional :: lower_tail, log_probability
      real(dp) :: high, log_lower, log_upper, low, midpoint
      logical :: invert_lower, lower, probability_is_log
      integer :: iteration

      lower = .true.
      if (present(lower_tail)) lower = lower_tail
      probability_is_log = .false.
      if (present(log_probability)) probability_is_log = log_probability
      if (degrees_freedom1 <= 0.0_dp .or. degrees_freedom2 <= 0.0_dp .or. &
          .not. ieee_is_finite(degrees_freedom1) .or. .not. ieee_is_finite(degrees_freedom2)) then
         value = ieee_value(0.0_dp,ieee_quiet_nan)
         return
      end if
      call probability_logs(probability,lower,probability_is_log,log_lower,log_upper)
      if (log_lower /= log_lower .or. log_upper /= log_upper) then
         value = ieee_value(0.0_dp,ieee_quiet_nan)
         return
      end if
      if (log_lower == ieee_value(0.0_dp,ieee_negative_inf)) then
         value = 0.0_dp
         return
      else if (log_upper == ieee_value(0.0_dp,ieee_negative_inf)) then
         value = ieee_value(0.0_dp,ieee_positive_inf)
         return
      end if

      invert_lower = log_lower <= log_upper
      low = 0.0_dp
      high = 1.0_dp
      if (invert_lower) then
         do while (r_pf(high,degrees_freedom1,degrees_freedom2,log_probability=.true.) < log_lower)
            if (high >= 0.25_dp*huge(1.0_dp)) exit
            high = 2.0_dp*high
         end do
      else
         do while (r_pf(high,degrees_freedom1,degrees_freedom2,lower_tail=.false., &
            log_probability=.true.) > log_upper)
            if (high >= 0.25_dp*huge(1.0_dp)) exit
            high = 2.0_dp*high
         end do
      end if
      do iteration = 1, 240
         midpoint = 0.5_dp*(low + high)
         if (invert_lower) then
            if (r_pf(midpoint,degrees_freedom1,degrees_freedom2,log_probability=.true.) < log_lower) then
               low = midpoint
            else
               high = midpoint
            end if
         else
            if (r_pf(midpoint,degrees_freedom1,degrees_freedom2,lower_tail=.false., &
               log_probability=.true.) > log_upper) then
               low = midpoint
            else
               high = midpoint
            end if
         end if
         if (high - low <= 8.0_dp*epsilon(1.0_dp)*max(1.0_dp,midpoint)) exit
      end do
      value = 0.5_dp*(low + high)
   end function r_qf

   pure elemental real(dp) function probability_result(lower_probability,upper_probability, &
      lower_tail,log_probability) result(value)
      real(dp), intent(in) :: lower_probability,upper_probability
      logical, intent(in) :: lower_tail,log_probability
      real(dp) :: selected

      selected = merge(lower_probability,upper_probability,lower_tail)
      if (log_probability) then
         if (lower_tail) then
            if (lower_probability == 0.0_dp) then
               value = ieee_value(0.0_dp,ieee_negative_inf)
            else if (upper_probability == 0.0_dp) then
               value = 0.0_dp
            else if (lower_probability > 0.5_dp) then
               value = r_log1p(-upper_probability)
            else
               value = log(lower_probability)
            end if
         else
            if (upper_probability == 0.0_dp) then
               value = ieee_value(0.0_dp,ieee_negative_inf)
            else if (lower_probability == 0.0_dp) then
               value = 0.0_dp
            else if (upper_probability > 0.5_dp) then
               value = r_log1p(-lower_probability)
            else
               value = log(upper_probability)
            end if
         end if
      else
         value = selected
      end if
   end function probability_result

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
