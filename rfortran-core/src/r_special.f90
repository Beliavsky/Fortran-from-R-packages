! SPDX-License-Identifier: MIT
module r_special
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_negative_inf, &
      ieee_quiet_nan, ieee_value
   use r_kinds, only : dp
   implicit none
   private

   public :: r_digamma, r_trigamma, r_log_beta, r_log_factorial, r_log_choose
   public :: r_regularized_beta, r_regularized_gamma_p, r_regularized_gamma_q

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

   pure elemental real(dp) function r_log_factorial(n) result(value)
      integer, intent(in) :: n

      if (n < 0) then
         value = ieee_value(0.0_dp, ieee_quiet_nan)
      else
         value = log_gamma(real(n + 1, dp))
      end if
   end function r_log_factorial

   pure elemental real(dp) function r_log_choose(n, k) result(value)
      integer, intent(in) :: n, k
      integer :: i, smaller

      if (n < 0 .or. k < 0 .or. k > n) then
         value = ieee_value(0.0_dp, ieee_negative_inf)
         return
      end if

      smaller = min(k, n - k)
      if (smaller <= 64) then
         value = 0.0_dp
         do i = 1, smaller
            value = value + log(real(n - smaller + i, dp)) - log(real(i, dp))
         end do
      else
         value = log_gamma(real(n + 1, dp)) - log_gamma(real(smaller + 1, dp)) - &
            log_gamma(real(n - smaller + 1, dp))
      end if
   end function r_log_choose

   pure elemental real(dp) function r_regularized_gamma_p(shape, x) result(value)
      real(dp), intent(in) :: shape, x

      if (.not. valid_gamma_arguments(shape, x)) then
         value = ieee_value(0.0_dp, ieee_quiet_nan)
      else if (x == 0.0_dp) then
         value = 0.0_dp
      else if (.not. ieee_is_finite(x)) then
         value = 1.0_dp
      else if (x < shape + 1.0_dp) then
         value = gamma_series(shape, x)
      else
         value = 1.0_dp - gamma_continued_fraction(shape, x)
      end if
      if (value == value) value = min(1.0_dp, max(0.0_dp, value))
   end function r_regularized_gamma_p

   pure elemental real(dp) function r_regularized_gamma_q(shape, x) result(value)
      real(dp), intent(in) :: shape, x

      if (.not. valid_gamma_arguments(shape, x)) then
         value = ieee_value(0.0_dp, ieee_quiet_nan)
      else if (x == 0.0_dp) then
         value = 1.0_dp
      else if (.not. ieee_is_finite(x)) then
         value = 0.0_dp
      else if (x < shape + 1.0_dp) then
         value = 1.0_dp - gamma_series(shape, x)
      else
         value = gamma_continued_fraction(shape, x)
      end if
      if (value == value) value = min(1.0_dp, max(0.0_dp, value))
   end function r_regularized_gamma_q

   pure elemental real(dp) function r_regularized_beta(x, shape1, shape2) result(value)
      real(dp), intent(in) :: x, shape1, shape2
      real(dp) :: factor

      if (x /= x .or. x < 0.0_dp .or. x > 1.0_dp) then
         value = ieee_value(0.0_dp, ieee_quiet_nan)
         return
      end if
      if (shape1 <= 0.0_dp .or. shape2 <= 0.0_dp) then
         value = ieee_value(0.0_dp, ieee_quiet_nan)
         return
      end if
      if (.not. ieee_is_finite(shape1) .or. .not. ieee_is_finite(shape2)) then
         value = ieee_value(0.0_dp, ieee_quiet_nan)
         return
      end if
      if (x == 0.0_dp) then
         value = 0.0_dp
         return
      end if
      if (x == 1.0_dp) then
         value = 1.0_dp
         return
      end if

      factor = exp(shape1*log(x) + shape2*log(1.0_dp - x) - r_log_beta(shape1, shape2))
      if (x < (shape1 + 1.0_dp)/(shape1 + shape2 + 2.0_dp)) then
         value = factor*beta_continued_fraction(shape1, shape2, x)/shape1
      else
         value = 1.0_dp - factor*beta_continued_fraction(shape2, shape1, 1.0_dp - x)/shape2
      end if
      value = min(1.0_dp, max(0.0_dp, value))
   end function r_regularized_beta

   pure elemental logical function valid_gamma_arguments(shape, x) result(valid)
      real(dp), intent(in) :: shape, x
      valid = shape > 0.0_dp .and. ieee_is_finite(shape) .and. x >= 0.0_dp .and. x == x
   end function valid_gamma_arguments

   pure elemental real(dp) function gamma_series(shape, x) result(value)
      real(dp), intent(in) :: shape, x
      real(dp) :: increment, shifted_shape, total
      integer :: iteration

      shifted_shape = shape
      total = 1.0_dp/shape
      increment = total
      do iteration = 1, 10000
         shifted_shape = shifted_shape + 1.0_dp
         increment = increment*x/shifted_shape
         total = total + increment
         if (abs(increment) <= 8.0_dp*epsilon(1.0_dp)*abs(total)) exit
      end do
      value = total*exp(-x + shape*log(x) - log_gamma(shape))
   end function gamma_series

   pure elemental real(dp) function gamma_continued_fraction(shape, x) result(value)
      real(dp), intent(in) :: shape, x
      real(dp), parameter :: floor_value = tiny(1.0_dp)/epsilon(1.0_dp)
      real(dp) :: coefficient, c, d, delta, denominator, h
      integer :: iteration

      denominator = x + 1.0_dp - shape
      if (abs(denominator) < floor_value) denominator = floor_value
      d = 1.0_dp/denominator
      c = 1.0_dp/floor_value
      h = d
      do iteration = 1, 10000
         coefficient = -real(iteration, dp)*(real(iteration, dp) - shape)
         denominator = denominator + 2.0_dp
         d = coefficient*d + denominator
         if (abs(d) < floor_value) d = sign(floor_value, d)
         c = denominator + coefficient/c
         if (abs(c) < floor_value) c = sign(floor_value, c)
         d = 1.0_dp/d
         delta = d*c
         h = h*delta
         if (abs(delta - 1.0_dp) <= 8.0_dp*epsilon(1.0_dp)) exit
      end do
      value = exp(-x + shape*log(x) - log_gamma(shape))*h
   end function gamma_continued_fraction

   pure elemental real(dp) function beta_continued_fraction(shape1, shape2, x) result(value)
      real(dp), intent(in) :: shape1, shape2, x
      real(dp), parameter :: floor_value = tiny(1.0_dp)/epsilon(1.0_dp)
      real(dp) :: c, d, delta, h, numerator, qab, qam, qap
      integer :: iteration, twice_iteration

      qab = shape1 + shape2
      qap = shape1 + 1.0_dp
      qam = shape1 - 1.0_dp
      c = 1.0_dp
      d = 1.0_dp - qab*x/qap
      if (abs(d) < floor_value) d = sign(floor_value, d)
      d = 1.0_dp/d
      h = d
      do iteration = 1, 10000
         twice_iteration = 2*iteration
         numerator = real(iteration, dp)*(shape2 - real(iteration, dp))*x
         numerator = numerator/((qam + real(twice_iteration, dp))*(shape1 + real(twice_iteration, dp)))
         d = 1.0_dp + numerator*d
         if (abs(d) < floor_value) d = sign(floor_value, d)
         c = 1.0_dp + numerator/c
         if (abs(c) < floor_value) c = sign(floor_value, c)
         d = 1.0_dp/d
         h = h*d*c
         numerator = -(shape1 + real(iteration, dp))*(qab + real(iteration, dp))*x
         numerator = numerator/((shape1 + real(twice_iteration, dp))*(qap + real(twice_iteration, dp)))
         d = 1.0_dp + numerator*d
         if (abs(d) < floor_value) d = sign(floor_value, d)
         c = 1.0_dp + numerator/c
         if (abs(c) < floor_value) c = sign(floor_value, c)
         d = 1.0_dp/d
         delta = d*c
         h = h*delta
         if (abs(delta - 1.0_dp) <= 8.0_dp*epsilon(1.0_dp)) exit
      end do
      value = h
   end function beta_continued_fraction

end module r_special
