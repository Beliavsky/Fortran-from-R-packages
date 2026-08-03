! SPDX-License-Identifier: GPL-2.0-or-later
module gnorm_special
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_positive_inf
   use gnorm_kinds, only : dp
   implicit none
   private
   public :: regularized_gamma_p, regularized_gamma_q
   public :: inverse_regularized_gamma_p, log1m, one_minus_exp

contains

   pure real(dp) function regularized_gamma_p(a, x) result(p)
      real(dp), intent(in) :: a, x
      real(dp) :: gln

      if (a <= 0.0_dp .or. x < 0.0_dp) then
         p = 0.0_dp
      else if (x <= tiny(1.0_dp)) then
         p = 0.0_dp
      else if (x < a + 1.0_dp) then
         gln = log_gamma(a)
         p = gamma_series(a, x, gln)
      else
         gln = log_gamma(a)
         p = 1.0_dp - gamma_continued_fraction(a, x, gln)
      end if
      p = min(1.0_dp, max(0.0_dp, p))
   end function regularized_gamma_p

   pure real(dp) function regularized_gamma_q(a, x) result(q)
      real(dp), intent(in) :: a, x
      real(dp) :: gln

      if (a <= 0.0_dp .or. x < 0.0_dp) then
         q = 0.0_dp
      else if (x <= tiny(1.0_dp)) then
         q = 1.0_dp
      else if (x < a + 1.0_dp) then
         gln = log_gamma(a)
         q = 1.0_dp - gamma_series(a, x, gln)
      else
         gln = log_gamma(a)
         q = gamma_continued_fraction(a, x, gln)
      end if
      q = min(1.0_dp, max(0.0_dp, q))
   end function regularized_gamma_q

   pure real(dp) function inverse_regularized_gamma_p(a, probability) result(x)
      real(dp), intent(in) :: a, probability
      real(dp) :: lower, upper, middle, p_mid
      integer :: i

      if (a <= 0.0_dp .or. probability < 0.0_dp .or. probability > 1.0_dp) then
         x = -1.0_dp
         return
      end if
      if (probability <= 0.0_dp) then
         x = 0.0_dp
         return
      end if
      if (probability >= 1.0_dp) then
         x = ieee_value(1.0_dp, ieee_positive_inf)
         return
      end if

      lower = 0.0_dp
      upper = max(1.0_dp, a)
      do i = 1, 1024
         if (regularized_gamma_p(a, upper) >= probability) exit
         if (upper >= huge(1.0_dp) / 2.0_dp) then
            x = ieee_value(1.0_dp, ieee_positive_inf)
            return
         end if
         upper = 2.0_dp * upper
      end do

      do i = 1, 140
         middle = 0.5_dp * (lower + upper)
         p_mid = regularized_gamma_p(a, middle)
         if (p_mid < probability) then
            lower = middle
         else
            upper = middle
         end if
         if (upper - lower <= 8.0_dp * epsilon(1.0_dp) * max(1.0_dp, middle)) exit
      end do
      x = 0.5_dp * (lower + upper)
   end function inverse_regularized_gamma_p

   pure real(dp) function log1m(x) result(value)
      real(dp), intent(in) :: x
      real(dp) :: term, sum_value
      integer :: k

      if (x <= -1.0_dp) then
         value = log(1.0_dp - x)
      else if (x >= 1.0_dp) then
         value = -ieee_value(1.0_dp, ieee_positive_inf)
      else if (abs(x) > 1.0e-4_dp) then
         value = log(1.0_dp - x)
      else
         term = x
         sum_value = 0.0_dp
         do k = 1, 12
            sum_value = sum_value - term / real(k, dp)
            term = term * x
         end do
         value = sum_value
      end if
   end function log1m

   pure real(dp) function one_minus_exp(x) result(value)
      real(dp), intent(in) :: x
      real(dp) :: term, sum_value
      integer :: k

      if (abs(x) > 1.0e-4_dp) then
         value = 1.0_dp - exp(x)
      else
         term = x
         sum_value = 0.0_dp
         do k = 1, 14
            sum_value = sum_value - term
            term = term * x / real(k + 1, dp)
         end do
         value = sum_value
      end if
   end function one_minus_exp

   pure real(dp) function gamma_series(a, x, gln) result(p)
      real(dp), intent(in) :: a, x, gln
      integer, parameter :: maxiter = 100000
      real(dp), parameter :: tolerance = 4.0e-15_dp
      real(dp) :: ap, delta, sum_value
      integer :: n

      ap = a
      sum_value = 1.0_dp / a
      delta = sum_value
      do n = 1, maxiter
         ap = ap + 1.0_dp
         delta = delta * x / ap
         sum_value = sum_value + delta
         if (abs(delta) <= abs(sum_value) * tolerance) exit
      end do
      p = sum_value * exp(-x + a * log(x) - gln)
   end function gamma_series

   pure real(dp) function gamma_continued_fraction(a, x, gln) result(q)
      real(dp), intent(in) :: a, x, gln
      integer, parameter :: maxiter = 100000
      real(dp), parameter :: tolerance = 4.0e-15_dp
      real(dp), parameter :: floor_value = tiny(1.0_dp) / tolerance
      real(dp) :: an, b, c, d, delta, h
      integer :: i

      b = x + 1.0_dp - a
      if (abs(b) < floor_value) b = sign(floor_value, b)
      c = 1.0_dp / floor_value
      d = 1.0_dp / b
      h = d
      do i = 1, maxiter
         an = -real(i, dp) * (real(i, dp) - a)
         b = b + 2.0_dp
         d = an * d + b
         if (abs(d) < floor_value) d = sign(floor_value, d)
         c = b + an / c
         if (abs(c) < floor_value) c = sign(floor_value, c)
         d = 1.0_dp / d
         delta = d * c
         h = h * delta
         if (abs(delta - 1.0_dp) <= tolerance) exit
      end do
      q = exp(-x + a * log(x) - gln) * h
   end function gamma_continued_fraction

end module gnorm_special
