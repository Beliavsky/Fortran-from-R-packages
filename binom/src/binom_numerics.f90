! SPDX-License-Identifier: GPL-2.0-or-later
! Numerical adapters for the binom translation; see NOTICE.md.
module binom_numerics
   use r_kinds, only: dp
   use r_special, only: r_regularized_beta, r_log_beta
   use r_transforms, only: r_log1p, r_expm1
   use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_is_finite, ieee_negative_inf
   implicit none
   private
   public :: dp, beta_probability, beta_quantile, log_binomial_mass, binomial_mass
   public :: binom_nan, valid_binomial, logistic, logit, binomial_deviance
contains
   pure elemental real(dp) function binom_nan() result(value)
      value = ieee_value(0.0_dp, ieee_quiet_nan)
   end function binom_nan

   pure elemental logical function valid_binomial(x, n, level) result(valid)
      real(dp), intent(in) :: x !! Number of successes, allowing fractional counts for power calculations.
      real(dp), intent(in) :: n !! Positive binomial trial count, at least as large as x.
      real(dp), intent(in) :: level !! Confidence probability strictly between zero and one.
      valid = .false.
      if (.not. ieee_is_finite(x) .or. .not. ieee_is_finite(n) .or. .not. ieee_is_finite(level)) return
      valid = n > 0.0_dp .and. x >= 0.0_dp .and. x <= n .and. level > 0.0_dp .and. level < 1.0_dp
   end function valid_binomial

   pure elemental real(dp) function beta_probability(x, a, b, upper_tail) result(value)
      real(dp), intent(in) :: x !! Beta variate, from zero through one.
      real(dp), intent(in) :: a !! Nonnegative first beta shape; zero uses the degenerate beta limit.
      real(dp), intent(in) :: b !! Nonnegative second beta shape; zero uses the degenerate beta limit.
      logical, intent(in), optional :: upper_tail !! If true compute the survival probability directly; default false.
      logical :: upper
      upper = .false.
      if (present(upper_tail)) upper = upper_tail
      if (a == 0.0_dp .or. b == 0.0_dp) then
         value = binom_nan()
         if (.not. all(ieee_is_finite([x, a, b]))) return
         if (x < 0.0_dp .or. x > 1.0_dp .or. a < 0.0_dp .or. b < 0.0_dp) return
         value = 0.0_dp
         if (x == 1.0_dp) then
            value = 1.0_dp
         else
            if (a == 0.0_dp) value = 1.0_dp
            if (a == 0.0_dp .and. b == 0.0_dp) value = 0.5_dp
         end if
         if (upper) value = 1.0_dp - value
         return
      end if
      if (upper) then
         value = r_regularized_beta(1.0_dp - x, b, a)
      else
         value = r_regularized_beta(x, a, b)
      end if
   end function beta_probability

   pure elemental real(dp) function beta_quantile(p, a, b) result(value)
      real(dp), intent(in) :: p !! Lower-tail beta probability in the closed unit interval.
      real(dp), intent(in) :: a !! Strictly positive first beta shape parameter.
      real(dp), intent(in) :: b !! Strictly positive second beta shape parameter.
      real(dp) :: target, aa, bb, low, high, mid, probability
      integer :: iteration
      value = binom_nan()
      if (.not. all(ieee_is_finite([p, a, b]))) return
      if (p < 0.0_dp .or. p > 1.0_dp .or. a <= 0.0_dp .or. b <= 0.0_dp) return
      if (p == 0.0_dp .or. p == 1.0_dp) then
         value = p
         return
      end if
      target = p
      aa = a
      bb = b
      if (p > 0.5_dp) then
         target = 1.0_dp - p
         aa = b
         bb = a
      end if
      low = 0.0_dp
      high = 1.0_dp
      do iteration = 1, 1200
         mid = low + 0.5_dp * (high - low)
         if (mid == low .or. mid == high) exit
         probability = beta_probability(mid, aa, bb)
         if (.not. ieee_is_finite(probability)) return
         if (probability < target) then
            low = mid
         else
            high = mid
         end if
         if (high - low <= 8.0_dp * epsilon(1.0_dp) * max(tiny(1.0_dp), mid)) exit
      end do
      value = low + 0.5_dp * (high - low)
      if (p > 0.5_dp) value = 1.0_dp - value
   end function beta_quantile

   pure elemental real(dp) function log_binomial_mass(x, n, p) result(value)
      real(dp), intent(in) :: x !! Success count in the closed interval zero through n.
      real(dp), intent(in) :: n !! Nonnegative trial count; fractional values use gamma interpolation.
      real(dp), intent(in) :: p !! Binomial success probability in the closed unit interval.
      value = binom_nan()
      if (.not. all(ieee_is_finite([x, n, p]))) return
      if (n < 0.0_dp .or. x < 0.0_dp .or. x > n .or. p < 0.0_dp .or. p > 1.0_dp) return
      value = ieee_value(0.0_dp, ieee_negative_inf)
      if (p == 0.0_dp) then
         if (x == 0.0_dp) value = 0.0_dp
      else if (p == 1.0_dp) then
         if (x == n) value = 0.0_dp
      else
         value = log_gamma(n + 1.0_dp) - log_gamma(x + 1.0_dp) - log_gamma(n - x + 1.0_dp) + &
            x * log(p) + (n - x) * r_log1p(-p)
      end if
   end function log_binomial_mass

   pure elemental real(dp) function binomial_mass(x, n, p) result(value)
      real(dp), intent(in) :: x !! Success count in the closed interval zero through n.
      real(dp), intent(in) :: n !! Nonnegative trial count.
      real(dp), intent(in) :: p !! Success probability in the closed unit interval.
      value = exp(log_binomial_mass(x, n, p))
   end function binomial_mass

   pure elemental real(dp) function logistic(x) result(value)
      real(dp), intent(in) :: x !! Real log odds; infinite values map to the corresponding probability endpoint.
      if (x >= 0.0_dp) then
         value = 1.0_dp / (1.0_dp + exp(-x))
      else
         value = exp(x) / (1.0_dp + exp(x))
      end if
   end function logistic

   pure elemental real(dp) function logit(p) result(value)
      real(dp), intent(in) :: p !! Probability strictly between zero and one.
      value = log(p) - r_log1p(-p)
   end function logit

   pure elemental real(dp) function binomial_deviance(p, x, n) result(value)
      real(dp), intent(in) :: p !! Candidate success probability strictly between zero and one.
      real(dp), intent(in) :: x !! Observed success count between zero and n.
      real(dp), intent(in) :: n !! Positive trial count.
      real(dp) :: estimate
      estimate = x / n
      value = 0.0_dp
      if (x > 0.0_dp) value = x * (log(estimate) - log(p))
      if (x < n) value = value + (n - x) * (r_log1p(-estimate) - r_log1p(-p))
      value = max(0.0_dp, 2.0_dp * value)
   end function binomial_deviance
end module binom_numerics
