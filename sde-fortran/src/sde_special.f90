! This file is part of sde-fortran, a translation/adaptation of the R package
! sde 2.0.21 by Stefano Maria Iacus.
! Original package code Copyright (C) 2006 S. M. Iacus.
! SPDX-License-Identifier: GPL-2.0-or-later
module sde_special
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
   use sde_kinds, only : dp, pi, sqrt_two, log_two_pi, tiny_dp
   use sde_interfaces, only : state_function
   implicit none
   private

   public :: normal_pdf
   public :: normal_logpdf
   public :: normal_cdf
   public :: normal_quantile
   public :: lognormal_pdf
   public :: lognormal_cdf
   public :: lognormal_quantile
   public :: gamma_pdf
   public :: gamma_logpdf
   public :: gamma_cdf
   public :: gamma_quantile
   public :: chi_square_pdf
   public :: chi_square_cdf
   public :: chi_square_quantile
   public :: regularized_gamma_p
   public :: safe_expm1
   public :: safe_log1p
   public :: expm1_over_x
   public :: log_sum_exp
   public :: integrate_adaptive
   public :: nan_dp

contains

   pure function nan_dp() result(value)
      real(dp) :: value
      value = ieee_value(0.0_dp, ieee_quiet_nan)
   end function nan_dp

   pure function safe_expm1(x) result(value)
      real(dp), intent(in) :: x
      real(dp) :: value
      real(dp) :: x2

      if (abs(x) < 1.0e-5_dp) then
         x2 = x*x
         value = x + 0.5_dp*x2 + x*x2/6.0_dp + x2*x2/24.0_dp + x2*x2*x/120.0_dp
      else
         value = exp(x) - 1.0_dp
      end if
   end function safe_expm1

   pure function safe_log1p(x) result(value)
      real(dp), intent(in) :: x
      real(dp) :: value
      real(dp) :: x2

      if (x <= -1.0_dp) then
         value = nan_dp()
      else if (abs(x) < 1.0e-5_dp) then
         x2 = x*x
         value = x - 0.5_dp*x2 + x*x2/3.0_dp - 0.25_dp*x2*x2 + 0.2_dp*x2*x2*x
      else
         value = log(1.0_dp + x)
      end if
   end function safe_log1p

   pure function expm1_over_x(x) result(value)
      real(dp), intent(in) :: x
      real(dp) :: value

      if (abs(x) < 1.0e-7_dp) then
         value = 1.0_dp + 0.5_dp*x + x*x/6.0_dp + x*x*x/24.0_dp
      else
         value = safe_expm1(x)/x
      end if
   end function expm1_over_x

   pure function normal_logpdf(x, mean, sd) result(value)
      real(dp), intent(in) :: x
      real(dp), intent(in), optional :: mean
      real(dp), intent(in), optional :: sd
      real(dp) :: value
      real(dp) :: location, scale, z

      location = 0.0_dp
      scale = 1.0_dp
      if (present(mean)) location = mean
      if (present(sd)) scale = sd
      if (scale <= 0.0_dp) then
         value = nan_dp()
         return
      end if
      z = (x-location)/scale
      value = -0.5_dp*(log_two_pi + z*z) - log(scale)
   end function normal_logpdf

   pure function normal_pdf(x, mean, sd) result(value)
      real(dp), intent(in) :: x
      real(dp), intent(in), optional :: mean
      real(dp), intent(in), optional :: sd
      real(dp) :: value

      value = exp(normal_logpdf(x, mean, sd))
   end function normal_pdf

   pure function normal_cdf(x, mean, sd, lower_tail) result(value)
      real(dp), intent(in) :: x
      real(dp), intent(in), optional :: mean
      real(dp), intent(in), optional :: sd
      logical, intent(in), optional :: lower_tail
      real(dp) :: value
      real(dp) :: location, scale, z
      logical :: lower

      location = 0.0_dp
      scale = 1.0_dp
      lower = .true.
      if (present(mean)) location = mean
      if (present(sd)) scale = sd
      if (present(lower_tail)) lower = lower_tail
      if (scale <= 0.0_dp) then
         value = nan_dp()
         return
      end if
      z = (x-location)/(scale*sqrt_two)
      if (lower) then
         value = 0.5_dp*erfc(-z)
      else
         value = 0.5_dp*erfc(z)
      end if
   end function normal_cdf

   pure function normal_quantile(p, mean, sd, lower_tail) result(value)
      real(dp), intent(in) :: p
      real(dp), intent(in), optional :: mean
      real(dp), intent(in), optional :: sd
      logical, intent(in), optional :: lower_tail
      real(dp) :: value
      real(dp) :: location, scale, prob, q, r
      logical :: lower
      real(dp), parameter :: a(6) = [ &
         -3.969683028665376e+01_dp, 2.209460984245205e+02_dp, &
         -2.759285104469687e+02_dp, 1.383577518672690e+02_dp, &
         -3.066479806614716e+01_dp, 2.506628277459239e+00_dp ]
      real(dp), parameter :: b(5) = [ &
         -5.447609879822406e+01_dp, 1.615858368580409e+02_dp, &
         -1.556989798598866e+02_dp, 6.680131188771972e+01_dp, &
         -1.328068155288572e+01_dp ]
      real(dp), parameter :: c(6) = [ &
         -7.784894002430293e-03_dp, -3.223964580411365e-01_dp, &
         -2.400758277161838e+00_dp, -2.549732539343734e+00_dp, &
          4.374664141464968e+00_dp, 2.938163982698783e+00_dp ]
      real(dp), parameter :: d(4) = [ &
          7.784695709041462e-03_dp, 3.224671290700398e-01_dp, &
          2.445134137142996e+00_dp, 3.754408661907416e+00_dp ]
      real(dp), parameter :: p_low = 0.02425_dp
      real(dp), parameter :: p_high = 1.0_dp-p_low

      location = 0.0_dp
      scale = 1.0_dp
      lower = .true.
      if (present(mean)) location = mean
      if (present(sd)) scale = sd
      if (present(lower_tail)) lower = lower_tail
      if (scale <= 0.0_dp .or. p < 0.0_dp .or. p > 1.0_dp) then
         value = nan_dp()
         return
      end if
      prob = p
      if (.not. lower) prob = 1.0_dp-p
      if (prob <= 0.0_dp) then
         value = -huge(1.0_dp)
         return
      else if (prob >= 1.0_dp) then
         value = huge(1.0_dp)
         return
      end if

      if (prob < p_low) then
         q = sqrt(-2.0_dp*log(prob))
         value = (((((c(1)*q+c(2))*q+c(3))*q+c(4))*q+c(5))*q+c(6))/ &
            ((((d(1)*q+d(2))*q+d(3))*q+d(4))*q+1.0_dp)
      else if (prob <= p_high) then
         q = prob-0.5_dp
         r = q*q
         value = (((((a(1)*r+a(2))*r+a(3))*r+a(4))*r+a(5))*r+a(6))*q/ &
            (((((b(1)*r+b(2))*r+b(3))*r+b(4))*r+b(5))*r+1.0_dp)
      else
         q = sqrt(-2.0_dp*log(1.0_dp-prob))
         value = -(((((c(1)*q+c(2))*q+c(3))*q+c(4))*q+c(5))*q+c(6))/ &
            ((((d(1)*q+d(2))*q+d(3))*q+d(4))*q+1.0_dp)
      end if

      ! One Halley refinement substantially improves the tail accuracy.
      q = normal_cdf(value)-prob
      value = value-q/normal_pdf(value)/(1.0_dp+0.5_dp*value*q/normal_pdf(value))
      value = location+scale*value
   end function normal_quantile

   pure function lognormal_pdf(x, meanlog, sdlog) result(value)
      real(dp), intent(in) :: x
      real(dp), intent(in) :: meanlog
      real(dp), intent(in) :: sdlog
      real(dp) :: value

      if (x <= 0.0_dp) then
         value = 0.0_dp
      else
         value = exp(normal_logpdf(log(x), meanlog, sdlog))/x
      end if
   end function lognormal_pdf

   pure function lognormal_cdf(x, meanlog, sdlog, lower_tail) result(value)
      real(dp), intent(in) :: x
      real(dp), intent(in) :: meanlog
      real(dp), intent(in) :: sdlog
      logical, intent(in), optional :: lower_tail
      real(dp) :: value
      logical :: lower

      lower = .true.
      if (present(lower_tail)) lower = lower_tail
      if (x <= 0.0_dp) then
         if (lower) then
            value = 0.0_dp
         else
            value = 1.0_dp
         end if
      else
         value = normal_cdf(log(x), meanlog, sdlog, lower)
      end if
   end function lognormal_cdf

   pure function lognormal_quantile(p, meanlog, sdlog, lower_tail) result(value)
      real(dp), intent(in) :: p
      real(dp), intent(in) :: meanlog
      real(dp), intent(in) :: sdlog
      logical, intent(in), optional :: lower_tail
      real(dp) :: value

      value = exp(normal_quantile(p, meanlog, sdlog, lower_tail))
   end function lognormal_quantile

   pure function regularized_gamma_p(a, x) result(value)
      real(dp), intent(in) :: a
      real(dp), intent(in) :: x
      real(dp) :: value
      real(dp) :: ap, del, summation, b, c, d, h, an, q
      integer :: n
      integer, parameter :: max_iter = 10000
      real(dp), parameter :: eps = 5.0e-15_dp
      real(dp), parameter :: fpmin = tiny_dp/eps

      if (a <= 0.0_dp .or. x < 0.0_dp) then
         value = nan_dp()
         return
      end if
      if (x <= 0.0_dp) then
         value = 0.0_dp
         return
      end if

      if (x < a+1.0_dp) then
         ap = a
         summation = 1.0_dp/a
         del = summation
         do n = 1, max_iter
            ap = ap+1.0_dp
            del = del*x/ap
            summation = summation+del
            if (abs(del) <= abs(summation)*eps) exit
         end do
         value = summation*exp(-x+a*log(x)-log_gamma(a))
      else
         b = x+1.0_dp-a
         c = 1.0_dp/fpmin
         d = 1.0_dp/b
         h = d
         do n = 1, max_iter
            an = -real(n, dp)*(real(n, dp)-a)
            b = b+2.0_dp
            d = an*d+b
            if (abs(d) < fpmin) d = fpmin
            c = b+an/c
            if (abs(c) < fpmin) c = fpmin
            d = 1.0_dp/d
            del = d*c
            h = h*del
            if (abs(del-1.0_dp) <= eps) exit
         end do
         q = exp(-x+a*log(x)-log_gamma(a))*h
         value = 1.0_dp-q
      end if
      value = max(0.0_dp, min(1.0_dp, value))
   end function regularized_gamma_p

   pure function gamma_logpdf(x, shape, scale) result(value)
      real(dp), intent(in) :: x
      real(dp), intent(in) :: shape
      real(dp), intent(in) :: scale
      real(dp) :: value

      if (shape <= 0.0_dp .or. scale <= 0.0_dp) then
         value = nan_dp()
      else if (x < 0.0_dp) then
         value = -huge(1.0_dp)
      else if (x <= 0.0_dp) then
         if (shape < 1.0_dp) then
            value = huge(1.0_dp)
         else if (shape <= 1.0_dp) then
            value = -log(scale)
         else
            value = -huge(1.0_dp)
         end if
      else
         value = (shape-1.0_dp)*log(x)-x/scale-log_gamma(shape)-shape*log(scale)
      end if
   end function gamma_logpdf

   pure function gamma_pdf(x, shape, scale) result(value)
      real(dp), intent(in) :: x
      real(dp), intent(in) :: shape
      real(dp), intent(in) :: scale
      real(dp) :: value
      real(dp) :: log_value

      log_value = gamma_logpdf(x, shape, scale)
      if (log_value >= log(huge(1.0_dp))) then
         value = huge(1.0_dp)
      else if (log_value <= log(tiny_dp)) then
         value = 0.0_dp
      else
         value = exp(log_value)
      end if
   end function gamma_pdf

   pure function gamma_cdf(x, shape, scale, lower_tail) result(value)
      real(dp), intent(in) :: x
      real(dp), intent(in) :: shape
      real(dp), intent(in) :: scale
      logical, intent(in), optional :: lower_tail
      real(dp) :: value
      logical :: lower

      lower = .true.
      if (present(lower_tail)) lower = lower_tail
      if (shape <= 0.0_dp .or. scale <= 0.0_dp) then
         value = nan_dp()
      else if (x <= 0.0_dp) then
         if (lower) then
            value = 0.0_dp
         else
            value = 1.0_dp
         end if
      else
         value = regularized_gamma_p(shape, x/scale)
         if (.not. lower) value = 1.0_dp-value
      end if
   end function gamma_cdf

   pure function gamma_quantile(p, shape, scale, lower_tail) result(value)
      real(dp), intent(in) :: p
      real(dp), intent(in) :: shape
      real(dp), intent(in) :: scale
      logical, intent(in), optional :: lower_tail
      real(dp) :: value
      real(dp) :: prob, lower_x, upper_x, mid
      logical :: lower
      integer :: iter

      lower = .true.
      if (present(lower_tail)) lower = lower_tail
      if (shape <= 0.0_dp .or. scale <= 0.0_dp .or. p < 0.0_dp .or. p > 1.0_dp) then
         value = nan_dp()
         return
      end if
      prob = p
      if (.not. lower) prob = 1.0_dp-p
      if (prob <= 0.0_dp) then
         value = 0.0_dp
         return
      else if (prob >= 1.0_dp) then
         value = huge(1.0_dp)
         return
      end if

      lower_x = 0.0_dp
      upper_x = max(scale, scale*(shape+8.0_dp*sqrt(shape)+8.0_dp))
      do while (gamma_cdf(upper_x, shape, scale) < prob)
         upper_x = 2.0_dp*upper_x
         if (upper_x >= huge(1.0_dp)/4.0_dp) exit
      end do
      do iter = 1, 200
         mid = 0.5_dp*(lower_x+upper_x)
         if (gamma_cdf(mid, shape, scale) < prob) then
            lower_x = mid
         else
            upper_x = mid
         end if
         if (upper_x-lower_x <= 2.0e-13_dp*max(1.0_dp, mid)) exit
      end do
      value = 0.5_dp*(lower_x+upper_x)
   end function gamma_quantile

   pure function chi_square_pdf(x, df) result(value)
      real(dp), intent(in) :: x
      real(dp), intent(in) :: df
      real(dp) :: value
      value = gamma_pdf(x, 0.5_dp*df, 2.0_dp)
   end function chi_square_pdf

   pure function chi_square_cdf(x, df, lower_tail) result(value)
      real(dp), intent(in) :: x
      real(dp), intent(in) :: df
      logical, intent(in), optional :: lower_tail
      real(dp) :: value
      value = gamma_cdf(x, 0.5_dp*df, 2.0_dp, lower_tail)
   end function chi_square_cdf

   pure function chi_square_quantile(p, df, lower_tail) result(value)
      real(dp), intent(in) :: p
      real(dp), intent(in) :: df
      logical, intent(in), optional :: lower_tail
      real(dp) :: value
      value = gamma_quantile(p, 0.5_dp*df, 2.0_dp, lower_tail)
   end function chi_square_quantile

   pure function log_sum_exp(values) result(value)
      real(dp), intent(in) :: values(:)
      real(dp) :: value
      real(dp) :: maximum

      if (size(values) == 0) then
         value = -huge(1.0_dp)
         return
      end if
      maximum = maxval(values)
      if (maximum <= -huge(1.0_dp)/2.0_dp) then
         value = maximum
      else
         value = maximum+log(sum(exp(values-maximum)))
      end if
   end function log_sum_exp

   function integrate_adaptive(f, a, b, theta, tolerance, max_depth) result(value)
      procedure(state_function) :: f
      real(dp), intent(in) :: a
      real(dp), intent(in) :: b
      real(dp), intent(in) :: theta(:)
      real(dp), intent(in), optional :: tolerance
      integer, intent(in), optional :: max_depth
      real(dp) :: value
      real(dp) :: tol, fa, fb, fm, whole
      integer :: depth

      tol = 1.0e-9_dp
      depth = 20
      if (present(tolerance)) tol = tolerance
      if (present(max_depth)) depth = max_depth
      if (abs(a-b) <= epsilon(1.0_dp)*max(1.0_dp, abs(a), abs(b))) then
         value = 0.0_dp
         return
      end if
      fa = f(a, theta)
      fb = f(b, theta)
      fm = f(0.5_dp*(a+b), theta)
      whole = (b-a)*(fa+4.0_dp*fm+fb)/6.0_dp
      value = adaptive_step(a, b, fa, fb, fm, whole, tol, depth)

   contains

      recursive function adaptive_step(left, right, f_left, f_right, f_mid, previous, eps, remaining) result(ans)
         real(dp), intent(in) :: left, right, f_left, f_right, f_mid, previous, eps
         integer, intent(in) :: remaining
         real(dp) :: ans
         real(dp) :: midpoint, left_mid, right_mid, f_left_mid, f_right_mid
         real(dp) :: left_area, right_area, total

         midpoint = 0.5_dp*(left+right)
         left_mid = 0.5_dp*(left+midpoint)
         right_mid = 0.5_dp*(midpoint+right)
         f_left_mid = f(left_mid, theta)
         f_right_mid = f(right_mid, theta)
         left_area = (midpoint-left)*(f_left+4.0_dp*f_left_mid+f_mid)/6.0_dp
         right_area = (right-midpoint)*(f_mid+4.0_dp*f_right_mid+f_right)/6.0_dp
         total = left_area+right_area
         if (remaining <= 0 .or. abs(total-previous) <= 15.0_dp*eps) then
            ans = total+(total-previous)/15.0_dp
         else
            ans = adaptive_step(left, midpoint, f_left, f_mid, f_left_mid, left_area, 0.5_dp*eps, remaining-1)+ &
               adaptive_step(midpoint, right, f_mid, f_right, f_right_mid, right_area, 0.5_dp*eps, remaining-1)
         end if
      end function adaptive_step

   end function integrate_adaptive

end module sde_special
