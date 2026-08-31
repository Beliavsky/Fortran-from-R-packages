! Copyright (C) 2011-2026 R. H. B. Christensen
! Modern Fortran translation, 2026. Distributed under GPL-2.0-or-later.
module ordinal_special
   use ordinal_kinds, only : dp
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_value, ieee_quiet_nan
   implicit none
   private
   real(dp), parameter :: pi = acos(-1.0_dp)
   public :: normal_cdf, normal_pdf, regularized_gamma_p
contains
   pure elemental real(dp) function normal_pdf(x) result(y)
      real(dp), intent(in) :: x !! Standard-normal variate at which to evaluate the density.
      y = exp(-0.5_dp*x*x)/sqrt(2.0_dp*pi)
   end function normal_pdf

   pure elemental real(dp) function normal_cdf(x) result(y)
      real(dp), intent(in) :: x !! Standard-normal variate at which to evaluate the lower-tail CDF.
      y = 0.5_dp*erfc(-x/sqrt(2.0_dp))
   end function normal_cdf

   pure real(dp) function regularized_gamma_p(a, x) result(p)
      real(dp), intent(in) :: a !! Positive gamma shape parameter.
      real(dp), intent(in) :: x !! Nonnegative argument of the lower regularized incomplete gamma function.
      integer, parameter :: itmax = 400
      real(dp), parameter :: eps = 5.0e-15_dp
      real(dp), parameter :: fpmin = tiny(1.0_dp)/eps
      real(dp) :: ap, del, sumv, b, c, d, h, an
      integer :: i
      if (ieee_is_nan(a) .or. ieee_is_nan(x) .or. a <= 0.0_dp .or. x < 0.0_dp) then
         p = ieee_value(0.0_dp, ieee_quiet_nan)
         return
      end if
      if (x <= 0.0_dp) then
         p = 0.0_dp
         return
      end if
      if (x < a + 1.0_dp) then
         ap = a
         sumv = 1.0_dp/a
         del = sumv
         do i = 1, itmax
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
         do i = 1, itmax
            an = -real(i, dp)*(real(i, dp) - a)
            b = b + 2.0_dp
            d = an*d + b
            if (abs(d) < fpmin) d = fpmin
            c = b + an/c
            if (abs(c) < fpmin) c = fpmin
            d = 1.0_dp/d
            del = d*c
            h = h*del
            if (abs(del - 1.0_dp) <= eps) exit
         end do
         p = 1.0_dp - exp(-x + a*log(x) - log_gamma(a))*h
      end if
      p = max(0.0_dp, min(1.0_dp, p))
   end function regularized_gamma_p
end module ordinal_special
