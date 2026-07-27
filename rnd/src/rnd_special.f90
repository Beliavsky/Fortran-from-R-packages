! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from RND 1.2, Copyright (C) 2017 Kam Hamidieh.
module rnd_special
   use rnd_kinds, only : dp, pi
   implicit none
   private
   public :: normal_pdf, normal_cdf, lognormal_pdf
   public :: beta_function, beta_pdf, regularized_beta
   public :: is_finite_real

contains

   elemental real(dp) function normal_pdf(x) result(value)
      real(dp), intent(in) :: x
      value = exp(-0.5_dp*x*x)/sqrt(2.0_dp*pi)
   end function normal_pdf

   elemental real(dp) function normal_cdf(x) result(value)
      real(dp), intent(in) :: x
      value = 0.5_dp*erfc(-x/sqrt(2.0_dp))
   end function normal_cdf

   elemental real(dp) function lognormal_pdf(x, meanlog, sdlog) result(value)
      real(dp), intent(in) :: x, meanlog, sdlog
      if (x <= 0.0_dp .or. sdlog <= 0.0_dp) then
         value = 0.0_dp
      else
         value = exp(-0.5_dp*((log(x)-meanlog)/sdlog)**2)/(x*sdlog*sqrt(2.0_dp*pi))
      end if
   end function lognormal_pdf

   elemental real(dp) function beta_function(a, b) result(value)
      real(dp), intent(in) :: a, b
      if (a <= 0.0_dp .or. b <= 0.0_dp) then
         value = huge(1.0_dp)
      else
         value = exp(log_gamma(a) + log_gamma(b) - log_gamma(a+b))
      end if
   end function beta_function

   elemental real(dp) function beta_pdf(x, a, b) result(value)
      real(dp), intent(in) :: x, a, b
      real(dp) :: log_value
      if (a <= 0.0_dp .or. b <= 0.0_dp .or. x <= 0.0_dp .or. x >= 1.0_dp) then
         value = 0.0_dp
      else
         log_value = (a-1.0_dp)*log(x) + (b-1.0_dp)*log(1.0_dp-x) &
            - log_gamma(a) - log_gamma(b) + log_gamma(a+b)
         value = exp(log_value)
      end if
   end function beta_pdf

   elemental real(dp) function regularized_beta(x, a, b) result(value)
      real(dp), intent(in) :: x, a, b
      real(dp) :: bt
      if (a <= 0.0_dp .or. b <= 0.0_dp) then
         value = 0.0_dp
      else if (x <= 0.0_dp) then
         value = 0.0_dp
      else if (x >= 1.0_dp) then
         value = 1.0_dp
      else
         bt = exp(log_gamma(a+b)-log_gamma(a)-log_gamma(b) &
            + a*log(x) + b*log(1.0_dp-x))
         if (x < (a+1.0_dp)/(a+b+2.0_dp)) then
            value = bt*beta_continued_fraction(a, b, x)/a
         else
            value = 1.0_dp - bt*beta_continued_fraction(b, a, 1.0_dp-x)/b
         end if
         value = max(0.0_dp, min(1.0_dp, value))
      end if
   end function regularized_beta

   elemental logical function is_finite_real(x) result(ok)
      use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
      real(dp), intent(in) :: x
      ok = ieee_is_finite(x)
   end function is_finite_real

   elemental real(dp) function beta_continued_fraction(a, b, x) result(value)
      real(dp), intent(in) :: a, b, x
      integer, parameter :: max_iter = 300
      real(dp), parameter :: eps = 3.0e-14_dp
      real(dp), parameter :: fpmin = tiny(1.0_dp)/eps
      integer :: m, m2
      real(dp) :: aa, c, d, del, h, qab, qam, qap

      qab = a+b
      qap = a+1.0_dp
      qam = a-1.0_dp
      c = 1.0_dp
      d = 1.0_dp-qab*x/qap
      if (abs(d) < fpmin) d = fpmin
      d = 1.0_dp/d
      h = d
      do m = 1, max_iter
         m2 = 2*m
         aa = real(m,dp)*(b-real(m,dp))*x/((qam+real(m2,dp))*(a+real(m2,dp)))
         d = 1.0_dp+aa*d
         if (abs(d) < fpmin) d = fpmin
         c = 1.0_dp+aa/c
         if (abs(c) < fpmin) c = fpmin
         d = 1.0_dp/d
         h = h*d*c
         aa = -(a+real(m,dp))*(qab+real(m,dp))*x &
            /((a+real(m2,dp))*(qap+real(m2,dp)))
         d = 1.0_dp+aa*d
         if (abs(d) < fpmin) d = fpmin
         c = 1.0_dp+aa/c
         if (abs(c) < fpmin) c = fpmin
         d = 1.0_dp/d
         del = d*c
         h = h*del
         if (abs(del-1.0_dp) <= eps) exit
      end do
      value = h
   end function beta_continued_fraction

end module rnd_special
