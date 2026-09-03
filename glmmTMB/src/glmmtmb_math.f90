! SPDX-License-Identifier: AGPL-3.0-only
! Derived from glmmTMB 1.1.14 computational sources; see NOTICE.md.
module glmmtmb_math
   use glmmtmb_kinds, only: dp
   use, intrinsic :: ieee_arithmetic, only: ieee_positive_inf, ieee_quiet_nan, ieee_value
   implicit none
   private
   real(dp), parameter :: log_two = log(2.0_dp)
   public :: expm1_safe, invlogit, log1mexp, log1p_safe, logaddexp, logsubexp
   public :: logcosh_safe, lambert_w0, log_bell_number
contains
   pure elemental real(dp) function log1p_safe(x) result(ans)
      real(dp), intent(in) :: x !! Increment in log(1+x), valid for x greater than -1.
      real(dp) :: term
      integer :: k
      if (x <= -1.0_dp) then
         ans = ieee_value(ans, ieee_quiet_nan)
      else if (abs(x) > 1.0e-4_dp) then
         ans = log(1.0_dp + x)
      else
         ans = 0.0_dp
         term = x
         do k = 1, 18
            if (mod(k, 2) == 1) then
               ans = ans + term / real(k, dp)
            else
               ans = ans - term / real(k, dp)
            end if
            term = term * x
         end do
      end if
   end function log1p_safe

   pure elemental real(dp) function expm1_safe(x) result(ans)
      real(dp), intent(in) :: x !! Exponent in exp(x)-1, evaluated accurately near zero.
      real(dp) :: term
      integer :: k
      if (abs(x) > 1.0e-4_dp) then
         ans = exp(x) - 1.0_dp
      else
         ans = 0.0_dp
         term = 1.0_dp
         do k = 1, 18
            term = term * x / real(k, dp)
            ans = ans + term
         end do
      end if
   end function expm1_safe

   pure elemental real(dp) function invlogit(x) result(ans)
      real(dp), intent(in) :: x !! Log-odds to transform to a probability in [0,1].
      if (x >= 0.0_dp) then
         ans = 1.0_dp / (1.0_dp + exp(-x))
      else
         ans = exp(x) / (1.0_dp + exp(x))
      end if
   end function invlogit

   pure elemental real(dp) function logaddexp(a, b) result(ans)
      real(dp), intent(in) :: a !! First log-scale quantity.
      real(dp), intent(in) :: b !! Second log-scale quantity.
      real(dp) :: m
      m = max(a, b)
      if (m == -ieee_value(m, ieee_positive_inf)) then
         ans = m
      else
         ans = m + log1p_safe(exp(min(a, b) - m))
      end if
   end function logaddexp

   pure elemental real(dp) function logsubexp(a, b) result(ans)
      real(dp), intent(in) :: a !! Logarithm of the minuend, must be at least b.
      real(dp), intent(in) :: b !! Logarithm of the subtrahend, must not exceed a.
      if (b > a) then
         ans = ieee_value(ans, ieee_quiet_nan)
      else if (b == a) then
         ans = -ieee_value(ans, ieee_positive_inf)
      else
         ans = a + log1p_safe(-exp(b - a))
      end if
   end function logsubexp

   pure elemental real(dp) function log1mexp(x) result(ans)
      real(dp), intent(in) :: x !! Nonpositive log probability x for log(1-exp(x)).
      if (x > 0.0_dp) then
         ans = ieee_value(ans, ieee_quiet_nan)
      else if (x < -log_two) then
         ans = log1p_safe(-exp(x))
      else
         ans = log(-expm1_safe(x))
      end if
   end function log1mexp

   pure elemental real(dp) function logcosh_safe(x) result(ans)
      real(dp), intent(in) :: x !! Real argument for a stable log(cosh(x)) evaluation.
      real(dp) :: ax
      ax = abs(x)
      ans = ax + log1p_safe(exp(-2.0_dp * ax)) - log_two
   end function logcosh_safe

   pure elemental real(dp) function lambert_w0(x) result(w)
      real(dp), intent(in) :: x !! Nonnegative argument for the principal Lambert W branch.
      real(dp) :: ew, f, denom, step
      integer :: iter
      if (x < 0.0_dp) then
         w = ieee_value(w, ieee_quiet_nan)
         return
      else if (x == 0.0_dp) then
         w = 0.0_dp
         return
      end if
      if (x < 1.0_dp) then
         w = x
      else
         w = log(x)
         if (x > 3.0_dp) w = w - log(w)
      end if
      do iter = 1, 40
         ew = exp(w)
         f = w * ew - x
         denom = ew * (w + 1.0_dp) - (w + 2.0_dp) * f / (2.0_dp * w + 2.0_dp)
         step = f / denom
         w = w - step
         if (abs(step) <= 8.0_dp * epsilon(1.0_dp) * max(1.0_dp, abs(w))) exit
      end do
   end function lambert_w0

   pure real(dp) function log_bell_number(n) result(ans)
      integer, intent(in) :: n !! Nonnegative Bell-number index, corresponding to an integer count.
      real(dp), allocatable :: row(:), next_row(:)
      integer :: i, j
      if (n < 0) then
         ans = ieee_value(ans, ieee_quiet_nan)
         return
      else if (n <= 1) then
         ans = 0.0_dp
         return
      end if
      allocate(row(n), next_row(n))
      row = -ieee_value(ans, ieee_positive_inf)
      next_row = row
      row(1) = 0.0_dp
      do i = 1, n - 1
         next_row(1) = row(i)
         do j = 2, i + 1
            next_row(j) = logaddexp(row(j - 1), next_row(j - 1))
         end do
         row = next_row
      end do
      ans = next_row(n)
   end function log_bell_number
end module glmmtmb_math
