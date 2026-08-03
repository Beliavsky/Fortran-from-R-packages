! SPDX-License-Identifier: GPL-2.0-or-later
module gnorm_distribution
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_value, &
      ieee_quiet_nan, ieee_positive_inf, ieee_negative_inf, ieee_is_nan
   use gnorm_kinds, only : dp, i8
   use gnorm_status, only : gnorm_success, gnorm_invalid_argument
   use gnorm_special, only : regularized_gamma_p, regularized_gamma_q, &
      inverse_regularized_gamma_p, log1m, one_minus_exp
   use gnorm_rng, only : gnorm_rng_state
   implicit none
   private
   public :: dgnorm, pgnorm, qgnorm, rgnorm, rgnorm_fill
   public :: gnorm_mean, gnorm_variance, gnorm_parameters_valid

contains

   elemental logical function gnorm_parameters_valid(alpha, beta) result(valid)
      real(dp), intent(in) :: alpha, beta
      valid = ieee_is_finite(alpha) .and. ieee_is_finite(beta) .and. &
         alpha > 0.0_dp .and. beta > 0.0_dp
   end function gnorm_parameters_valid

   elemental real(dp) function dgnorm(x, mu, alpha, beta, log_density) result(value)
      real(dp), intent(in) :: x
      real(dp), intent(in), optional :: mu, alpha, beta
      logical, intent(in), optional :: log_density
      real(dp) :: location, scale, shape, log_value, power_value
      logical :: return_log

      location = 0.0_dp
      scale = 1.0_dp
      shape = 1.0_dp
      return_log = .false.
      if (present(mu)) location = mu
      if (present(alpha)) scale = alpha
      if (present(beta)) shape = beta
      if (present(log_density)) return_log = log_density

      if (.not. gnorm_parameters_valid(scale, shape) .or. .not. ieee_is_finite(location)) then
         value = ieee_value(1.0_dp, ieee_quiet_nan)
         return
      end if
      if (.not. ieee_is_finite(x)) then
         if (return_log) then
            value = ieee_value(1.0_dp, ieee_negative_inf)
         else
            value = 0.0_dp
         end if
         return
      end if

      power_value = standardized_power(abs(x - location), scale, shape)
      log_value = log(shape) - log(2.0_dp) - log(scale) - &
         log_gamma(1.0_dp / shape) - power_value
      if (return_log) then
         value = log_value
      else if (log_value < log(tiny(1.0_dp))) then
         value = 0.0_dp
      else
         value = exp(log_value)
      end if
   end function dgnorm

   elemental real(dp) function pgnorm(q, mu, alpha, beta, lower_tail, &
      log_probability) result(value)
      real(dp), intent(in) :: q
      real(dp), intent(in), optional :: mu, alpha, beta
      logical, intent(in), optional :: lower_tail, log_probability
      real(dp) :: location, scale, shape, t, tail_probability, probability
      logical :: lower, return_log, requested_is_small_tail

      location = 0.0_dp
      scale = 1.0_dp
      shape = 1.0_dp
      lower = .true.
      return_log = .false.
      if (present(mu)) location = mu
      if (present(alpha)) scale = alpha
      if (present(beta)) shape = beta
      if (present(lower_tail)) lower = lower_tail
      if (present(log_probability)) return_log = log_probability

      if (.not. gnorm_parameters_valid(scale, shape) .or. .not. ieee_is_finite(location)) then
         value = ieee_value(1.0_dp, ieee_quiet_nan)
         return
      end if
      if (.not. ieee_is_finite(q)) then
         if (ieee_is_nan(q)) then
            value = ieee_value(1.0_dp, ieee_quiet_nan)
         else if (q < 0.0_dp) then
            probability = merge(0.0_dp, 1.0_dp, lower)
            value = probability_result(probability, return_log)
         else
            probability = merge(1.0_dp, 0.0_dp, lower)
            value = probability_result(probability, return_log)
         end if
         return
      end if

      t = standardized_power(abs(q - location), scale, shape)
      if (ieee_is_finite(t)) then
         tail_probability = 0.5_dp * regularized_gamma_q(1.0_dp / shape, t)
      else
         tail_probability = 0.0_dp
      end if
      requested_is_small_tail = (q < location .and. lower) .or. &
         (q >= location .and. .not. lower)
      if (requested_is_small_tail) then
         probability = tail_probability
         if (return_log) then
            if (probability <= 0.0_dp) then
               value = ieee_value(1.0_dp, ieee_negative_inf)
            else
               value = log(probability)
            end if
         else
            value = probability
         end if
      else
         probability = 1.0_dp - tail_probability
         if (return_log) then
            value = log1m(tail_probability)
         else
            value = probability
         end if
      end if
   end function pgnorm

   elemental real(dp) function qgnorm(p, mu, alpha, beta, lower_tail, &
      log_probability) result(value)
      real(dp), intent(in) :: p
      real(dp), intent(in), optional :: mu, alpha, beta
      logical, intent(in), optional :: lower_tail, log_probability
      real(dp) :: location, scale, shape, lower_probability, y, gamma_quantile
      real(dp) :: deviation
      logical :: lower, input_is_log

      location = 0.0_dp
      scale = 1.0_dp
      shape = 1.0_dp
      lower = .true.
      input_is_log = .false.
      if (present(mu)) location = mu
      if (present(alpha)) scale = alpha
      if (present(beta)) shape = beta
      if (present(lower_tail)) lower = lower_tail
      if (present(log_probability)) input_is_log = log_probability

      if (.not. gnorm_parameters_valid(scale, shape) .or. .not. ieee_is_finite(location)) then
         value = ieee_value(1.0_dp, ieee_quiet_nan)
         return
      end if

      if (input_is_log) then
         if (p > 0.0_dp .or. ieee_is_nan(p)) then
            value = ieee_value(1.0_dp, ieee_quiet_nan)
            return
         end if
         if (lower) then
            lower_probability = exp(p)
         else
            lower_probability = one_minus_exp(p)
         end if
      else
         if (p < 0.0_dp .or. p > 1.0_dp .or. ieee_is_nan(p)) then
            value = ieee_value(1.0_dp, ieee_quiet_nan)
            return
         end if
         lower_probability = merge(p, 1.0_dp - p, lower)
      end if

      if (lower_probability <= 0.0_dp) then
         value = ieee_value(1.0_dp, ieee_negative_inf)
         return
      else if (lower_probability >= 1.0_dp) then
         value = ieee_value(1.0_dp, ieee_positive_inf)
         return
      else if (abs(lower_probability - 0.5_dp) <= epsilon(1.0_dp)) then
         value = location
         return
      end if

      y = 2.0_dp * abs(lower_probability - 0.5_dp)
      gamma_quantile = inverse_regularized_gamma_p(1.0_dp / shape, y)
      if (.not. ieee_is_finite(gamma_quantile)) then
         deviation = ieee_value(1.0_dp, ieee_positive_inf)
      else if (gamma_quantile <= tiny(1.0_dp)) then
         deviation = 0.0_dp
      else
         deviation = scale * exp(log(gamma_quantile) / shape)
      end if
      value = location + sign(deviation, lower_probability - 0.5_dp)
   end function qgnorm

   function rgnorm(n, mu, alpha, beta, seed, status) result(values)
      integer, intent(in) :: n
      real(dp), intent(in), optional :: mu, alpha, beta
      integer(i8), intent(in), optional :: seed
      integer, intent(out), optional :: status
      real(dp), allocatable :: values(:)

      allocate(values(max(0, n)))
      call rgnorm_fill(values, mu, alpha, beta, seed, status)
   end function rgnorm

   subroutine rgnorm_fill(values, mu, alpha, beta, seed, status)
      real(dp), intent(out) :: values(:)
      real(dp), intent(in), optional :: mu, alpha, beta
      integer(i8), intent(in), optional :: seed
      integer, intent(out), optional :: status
      type(gnorm_rng_state) :: rng
      real(dp) :: location, scale, shape, magnitude
      integer(i8) :: clock_count
      integer :: i

      location = 0.0_dp
      scale = 1.0_dp
      shape = 1.0_dp
      if (present(mu)) location = mu
      if (present(alpha)) scale = alpha
      if (present(beta)) shape = beta
      if (present(status)) status = gnorm_success

      if (.not. gnorm_parameters_valid(scale, shape) .or. .not. ieee_is_finite(location)) then
         values = ieee_value(1.0_dp, ieee_quiet_nan)
         if (present(status)) status = gnorm_invalid_argument
         return
      end if
      if (present(seed)) then
         call rng%seed(seed)
      else
         call system_clock(count=clock_count)
         call rng%seed(clock_count)
      end if

      do i = 1, size(values)
         magnitude = scale * rng%gamma(1.0_dp / shape)**(1.0_dp / shape)
         if (rng%uniform() < 0.5_dp) magnitude = -magnitude
         values(i) = location + magnitude
      end do
   end subroutine rgnorm_fill

   elemental real(dp) function gnorm_mean(mu) result(value)
      real(dp), intent(in), optional :: mu
      value = 0.0_dp
      if (present(mu)) value = mu
   end function gnorm_mean

   elemental real(dp) function gnorm_variance(alpha, beta) result(value)
      real(dp), intent(in), optional :: alpha, beta
      real(dp) :: scale, shape

      scale = 1.0_dp
      shape = 1.0_dp
      if (present(alpha)) scale = alpha
      if (present(beta)) shape = beta
      if (.not. gnorm_parameters_valid(scale, shape)) then
         value = ieee_value(1.0_dp, ieee_quiet_nan)
      else
         value = scale*scale * exp(log_gamma(3.0_dp / shape) - &
            log_gamma(1.0_dp / shape))
      end if
   end function gnorm_variance

   elemental real(dp) function standardized_power(distance, scale, shape) result(value)
      real(dp), intent(in) :: distance, scale, shape
      real(dp) :: log_value

      if (distance <= tiny(1.0_dp)) then
         value = 0.0_dp
         return
      end if
      log_value = shape * (log(distance) - log(scale))
      if (log_value >= log(huge(1.0_dp))) then
         value = ieee_value(1.0_dp, ieee_positive_inf)
      else if (log_value <= log(tiny(1.0_dp))) then
         value = 0.0_dp
      else
         value = exp(log_value)
      end if
   end function standardized_power

   elemental real(dp) function probability_result(probability, return_log) result(value)
      real(dp), intent(in) :: probability
      logical, intent(in) :: return_log

      if (return_log) then
         if (probability <= 0.0_dp) then
            value = ieee_value(1.0_dp, ieee_negative_inf)
         else
            value = log(probability)
         end if
      else
         value = probability
      end if
   end function probability_result

end module gnorm_distribution
