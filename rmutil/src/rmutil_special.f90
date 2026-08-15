! rmutil computational translation
! Copyright (C) 1998-2001 J.K. Lindsey
! Copyright (C) 2026 OpenAI (modern Fortran translation)
! SPDX-License-Identifier: GPL-2.0-or-later
module rmutil_special
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, ieee_positive_inf
   use rmutil_kinds, only : dp, pi, sqrt2, log2pi
   implicit none
   private
   public :: normal_pdf, normal_cdf, normal_quantile
   public :: regularized_gamma_p, regularized_gamma_q, gamma_quantile
   public :: log_beta, log_choose, poisson_pmf, poisson_cdf
   public :: negative_binomial_pmf, negative_binomial_cdf
   public :: random_normal, random_gamma, random_integer_quantile
   public :: log1p_r, expm1_r, bisection_quantile

   abstract interface
      function cdf_scalar(x) result(p)
         import dp
         real(dp), intent(in) :: x
         real(dp) :: p
      end function cdf_scalar
      function discrete_cdf(k) result(p)
         import dp
         integer, intent(in) :: k
         real(dp) :: p
      end function discrete_cdf
   end interface

contains

   elemental real(dp) function nan_dp() result(x)
      x = ieee_value(0.0_dp, ieee_quiet_nan)
   end function nan_dp

   elemental real(dp) function inf_dp() result(x)
      x = ieee_value(0.0_dp, ieee_positive_inf)
   end function inf_dp

   elemental real(dp) function log1p_r(x) result(y)
      real(dp), intent(in) :: x
      if (abs(x) < 1.0e-8_dp) then
         y = x - x*x/2.0_dp + x**3/3.0_dp - x**4/4.0_dp
      else
         y = log(1.0_dp + x)
      end if
   end function log1p_r

   elemental real(dp) function expm1_r(x) result(y)
      real(dp), intent(in) :: x
      if (abs(x) < 1.0e-6_dp) then
         y = x + x*x/2.0_dp + x**3/6.0_dp + x**4/24.0_dp
      else
         y = exp(x) - 1.0_dp
      end if
   end function expm1_r

   elemental real(dp) function normal_pdf(x) result(y)
      real(dp), intent(in) :: x
      y = exp(-0.5_dp*x*x - 0.5_dp*log2pi)
   end function normal_pdf

   elemental real(dp) function normal_cdf(x) result(p)
      real(dp), intent(in) :: x
      p = 0.5_dp*erfc(-x/sqrt2)
   end function normal_cdf

   elemental real(dp) function normal_quantile(p) result(x)
      real(dp), intent(in) :: p
      real(dp), parameter :: a(6) = [ &
         -3.969683028665376e1_dp, 2.209460984245205e2_dp, &
         -2.759285104469687e2_dp, 1.383577518672690e2_dp, &
         -3.066479806614716e1_dp, 2.506628277459239_dp ]
      real(dp), parameter :: b(5) = [ &
         -5.447609879822406e1_dp, 1.615858368580409e2_dp, &
         -1.556989798598866e2_dp, 6.680131188771972e1_dp, &
         -1.328068155288572e1_dp ]
      real(dp), parameter :: c(6) = [ &
         -7.784894002430293e-3_dp, -3.223964580411365e-1_dp, &
         -2.400758277161838_dp, -2.549732539343734_dp, &
          4.374664141464968_dp, 2.938163982698783_dp ]
      real(dp), parameter :: d(4) = [ &
          7.784695709041462e-3_dp, 3.224671290700398e-1_dp, &
          2.445134137142996_dp, 3.754408661907416_dp ]
      real(dp), parameter :: plow = 0.02425_dp, phigh = 1.0_dp - plow
      real(dp) :: q, r, e, u
      if (p <= 0.0_dp) then
         x = merge(-inf_dp(), nan_dp(), p == 0.0_dp)
         return
      else if (p >= 1.0_dp) then
         x = merge(inf_dp(), nan_dp(), p == 1.0_dp)
         return
      end if
      if (p < plow) then
         q = sqrt(-2.0_dp*log(p))
         x = (((((c(1)*q+c(2))*q+c(3))*q+c(4))*q+c(5))*q+c(6)) / &
            ((((d(1)*q+d(2))*q+d(3))*q+d(4))*q+1.0_dp)
      else if (p <= phigh) then
         q = p - 0.5_dp
         r = q*q
         x = (((((a(1)*r+a(2))*r+a(3))*r+a(4))*r+a(5))*r+a(6))*q / &
            (((((b(1)*r+b(2))*r+b(3))*r+b(4))*r+b(5))*r+1.0_dp)
      else
         q = sqrt(-2.0_dp*log(1.0_dp-p))
         x = -(((((c(1)*q+c(2))*q+c(3))*q+c(4))*q+c(5))*q+c(6)) / &
            ((((d(1)*q+d(2))*q+d(3))*q+d(4))*q+1.0_dp)
      end if
      e = normal_cdf(x) - p
      u = e/max(normal_pdf(x), tiny(1.0_dp))
      x = x - u/(1.0_dp + 0.5_dp*x*u)
   end function normal_quantile

   pure real(dp) function regularized_gamma_p(a, x) result(p)
      real(dp), intent(in) :: a, x
      integer, parameter :: itmax = 10000
      real(dp), parameter :: eps = 4.0e-15_dp
      real(dp), parameter :: fpmin = tiny(1.0_dp)/eps
      integer :: n
      real(dp) :: ap, del, sumv, b, c, d, h, an, q
      if (a <= 0.0_dp .or. x < 0.0_dp) then
         p = nan_dp()
         return
      else if (x == 0.0_dp) then
         p = 0.0_dp
         return
      end if
      if (x < a + 1.0_dp) then
         ap = a
         sumv = 1.0_dp/a
         del = sumv
         do n = 1, itmax
            ap = ap + 1.0_dp
            del = del*x/ap
            sumv = sumv + del
            if (abs(del) <= abs(sumv)*eps) exit
         end do
         p = sumv*exp(-x + a*log(x) - log_gamma(a))
      else
         b = x + 1.0_dp - a
         c = 1.0_dp/fpmin
         d = 1.0_dp/b
         h = d
         do n = 1, itmax
            an = -real(n,dp)*(real(n,dp)-a)
            b = b + 2.0_dp
            d = an*d + b
            if (abs(d) < fpmin) d = fpmin
            c = b + an/c
            if (abs(c) < fpmin) c = fpmin
            d = 1.0_dp/d
            del = d*c
            h = h*del
            if (abs(del-1.0_dp) <= eps) exit
         end do
         q = exp(-x + a*log(x) - log_gamma(a))*h
         p = max(0.0_dp, min(1.0_dp, 1.0_dp-q))
      end if
   end function regularized_gamma_p

   pure real(dp) function regularized_gamma_q(a, x) result(q)
      real(dp), intent(in) :: a, x
      q = 1.0_dp - regularized_gamma_p(a, x)
   end function regularized_gamma_q

   pure real(dp) function gamma_quantile(p, shape, scale) result(x)
      real(dp), intent(in) :: p, shape, scale
      real(dp) :: lo, hi, mid, cdf
      integer :: iter
      if (p <= 0.0_dp) then
         x = 0.0_dp
         return
      else if (p >= 1.0_dp) then
         x = inf_dp()
         return
      end if
      lo = 0.0_dp
      hi = max(scale*shape, scale)
      do while (regularized_gamma_p(shape, hi/scale) < p)
         hi = 2.0_dp*hi
         if (hi > huge(1.0_dp)/4.0_dp) exit
      end do
      do iter = 1, 160
         mid = 0.5_dp*(lo+hi)
         cdf = regularized_gamma_p(shape, mid/scale)
         if (cdf < p) then
            lo = mid
         else
            hi = mid
         end if
         if (hi-lo <= 2.0e-13_dp*max(1.0_dp,mid)) exit
      end do
      x = 0.5_dp*(lo+hi)
   end function gamma_quantile

   elemental real(dp) function log_beta(a, b) result(v)
      real(dp), intent(in) :: a, b
      v = log_gamma(a) + log_gamma(b) - log_gamma(a+b)
   end function log_beta

   elemental real(dp) function log_choose(n, k) result(v)
      integer, intent(in) :: n, k
      if (k < 0 .or. k > n .or. n < 0) then
         v = -inf_dp()
      else
         v = log_gamma(real(n+1,dp)) - log_gamma(real(k+1,dp)) - &
            log_gamma(real(n-k+1,dp))
      end if
   end function log_choose

   elemental real(dp) function poisson_pmf(k, mu) result(p)
      integer, intent(in) :: k
      real(dp), intent(in) :: mu
      if (k < 0 .or. mu < 0.0_dp) then
         p = 0.0_dp
      else if (mu == 0.0_dp) then
         p = merge(1.0_dp, 0.0_dp, k == 0)
      else
         p = exp(real(k,dp)*log(mu) - mu - log_gamma(real(k+1,dp)))
      end if
   end function poisson_pmf

   real(dp) function poisson_cdf(k, mu) result(p)
      integer, intent(in) :: k
      real(dp), intent(in) :: mu
      if (k < 0) then
         p = 0.0_dp
      else if (mu == 0.0_dp) then
         p = 1.0_dp
      else
         p = regularized_gamma_q(real(k+1,dp), mu)
      end if
   end function poisson_cdf

   elemental real(dp) function negative_binomial_pmf(k, size, prob) result(p)
      integer, intent(in) :: k
      real(dp), intent(in) :: size, prob
      if (k < 0 .or. size <= 0.0_dp .or. prob <= 0.0_dp .or. prob > 1.0_dp) then
         p = 0.0_dp
      else
         p = exp(log_gamma(real(k,dp)+size) - log_gamma(size) - &
            log_gamma(real(k+1,dp)) + size*log(prob) + real(k,dp)*log(1.0_dp-prob))
      end if
   end function negative_binomial_pmf

   real(dp) function negative_binomial_cdf(k, size, prob) result(p)
      integer, intent(in) :: k
      real(dp), intent(in) :: size, prob
      integer :: j
      if (k < 0) then
         p = 0.0_dp
      else
         p = 0.0_dp
         do j = 0, k
            p = p + negative_binomial_pmf(j, size, prob)
         end do
         p = min(1.0_dp, p)
      end if
   end function negative_binomial_cdf

   real(dp) function bisection_quantile(p, f, lo0, hi0) result(x)
      real(dp), intent(in) :: p, lo0, hi0
      procedure(cdf_scalar) :: f
      real(dp) :: lo, hi, mid
      integer :: iter
      if (p <= 0.0_dp) then
         x = lo0
         return
      else if (p >= 1.0_dp) then
         x = hi0
         return
      end if
      lo = lo0
      hi = hi0
      do iter = 1, 180
         mid = 0.5_dp*(lo+hi)
         if (f(mid) < p) then
            lo = mid
         else
            hi = mid
         end if
         if (hi-lo <= 2.0e-13_dp*max(1.0_dp,abs(mid))) exit
      end do
      x = 0.5_dp*(lo+hi)
   end function bisection_quantile

   integer function random_integer_quantile(f) result(k)
      procedure(discrete_cdf) :: f
      real(dp) :: u
      integer :: hi, lo, mid
      call random_number(u)
      lo = 0
      hi = 16
      do
         if (f(hi) >= u) exit
         if (hi >= huge(hi)/4) exit
         hi = 2*hi
      end do
      do while (lo < hi)
         mid = (lo+hi)/2
         if (f(mid) >= u) then
            hi = mid
         else
            lo = mid + 1
         end if
      end do
      k = lo
   end function random_integer_quantile

   real(dp) function random_normal() result(z)
      real(dp) :: u1, u2
      call random_number(u1)
      call random_number(u2)
      u1 = max(u1, tiny(1.0_dp))
      z = sqrt(-2.0_dp*log(u1))*cos(2.0_dp*pi*u2)
   end function random_normal

   recursive real(dp) function random_gamma(shape, scale) result(x)
      real(dp), intent(in) :: shape, scale
      real(dp) :: d, c, z, u, v
      if (shape <= 0.0_dp .or. scale <= 0.0_dp) then
         x = nan_dp()
         return
      end if
      if (shape < 1.0_dp) then
         call random_number(u)
         x = random_gamma(shape+1.0_dp, scale)*u**(1.0_dp/shape)
         return
      end if
      d = shape - 1.0_dp/3.0_dp
      c = 1.0_dp/sqrt(9.0_dp*d)
      do
         do
            z = random_normal()
            v = 1.0_dp + c*z
            if (v > 0.0_dp) exit
         end do
         v = v**3
         call random_number(u)
         if (u < 1.0_dp - 0.0331_dp*z**4) exit
         if (log(u) < 0.5_dp*z*z + d*(1.0_dp-v+log(v))) exit
      end do
      x = scale*d*v
   end function random_gamma

end module rmutil_special
