! SPDX-License-Identifier: GPL-2.0-only
module poibin_special
   use poibin_kinds, only : dp
   implicit none
   private
   public :: normal_cdf, normal_pdf, poisson_cdf

contains

   pure real(dp) function normal_pdf(x) result(y)
      real(dp), intent(in) :: x
      real(dp), parameter :: inv_sqrt_2pi = 0.398942280401432677939946059934_dp
      y = inv_sqrt_2pi * exp(-0.5_dp*x*x)
   end function normal_pdf

   pure real(dp) function normal_cdf(x) result(y)
      real(dp), intent(in) :: x
      y = 0.5_dp * erfc(-x / sqrt(2.0_dp))
   end function normal_cdf

   pure real(dp) function poisson_cdf(k, lambda) result(y)
      integer, intent(in) :: k
      real(dp), intent(in) :: lambda

      if (k < 0) then
         y = 0.0_dp
      else if (lambda < 0.0_dp) then
         y = 0.0_dp
      else if (lambda <= 0.0_dp) then
         y = 1.0_dp
      else
         y = regularized_gamma_q(real(k + 1, dp), lambda)
         y = max(0.0_dp, min(1.0_dp, y))
      end if
   end function poisson_cdf

   pure real(dp) function regularized_gamma_q(a, x) result(q)
      real(dp), intent(in) :: a, x
      real(dp) :: p

      if (a <= 0.0_dp .or. x < 0.0_dp) then
         q = 0.0_dp
      else if (x <= 0.0_dp) then
         q = 1.0_dp
      else if (x < a + 1.0_dp) then
         p = gamma_series(a, x)
         q = 1.0_dp - p
      else
         q = gamma_contfrac(a, x)
      end if
      q = max(0.0_dp, min(1.0_dp, q))
   end function regularized_gamma_q

   pure real(dp) function gamma_series(a, x) result(p)
      real(dp), intent(in) :: a, x
      integer, parameter :: maxit = 10000
      real(dp), parameter :: eps = 8.0_dp*epsilon(1.0_dp)
      real(dp) :: ap, del, sumv
      integer :: n

      if (x <= 0.0_dp) then
         p = 0.0_dp
         return
      end if

      ap = a
      sumv = 1.0_dp/a
      del = sumv
      do n = 1, maxit
         ap = ap + 1.0_dp
         del = del*x/ap
         sumv = sumv + del
         if (abs(del) <= abs(sumv)*eps) exit
      end do
      p = sumv*exp(-x + a*log(x) - log_gamma(a))
   end function gamma_series

   pure real(dp) function gamma_contfrac(a, x) result(q)
      real(dp), intent(in) :: a, x
      integer, parameter :: maxit = 10000
      real(dp), parameter :: eps = 8.0_dp*epsilon(1.0_dp)
      real(dp), parameter :: fpmin = tiny(1.0_dp)/epsilon(1.0_dp)
      real(dp) :: an, b, c, d, del, h
      integer :: i

      b = x + 1.0_dp - a
      c = 1.0_dp/fpmin
      d = 1.0_dp/max(abs(b), fpmin)
      if (b < 0.0_dp) d = -d
      h = d
      do i = 1, maxit
         an = -real(i, dp)*(real(i, dp) - a)
         b = b + 2.0_dp
         d = an*d + b
         if (abs(d) < fpmin) d = sign(fpmin, d)
         c = b + an/c
         if (abs(c) < fpmin) c = sign(fpmin, c)
         d = 1.0_dp/d
         del = d*c
         h = h*del
         if (abs(del - 1.0_dp) <= eps) exit
      end do
      q = exp(-x + a*log(x) - log_gamma(a))*h
   end function gamma_contfrac

end module poibin_special
