! Modern Fortran translation of R package skewunit.
! SPDX-License-Identifier: GPL-2.0-or-later
module skewunit_special
   use skewunit_kinds, only : dp, sqrt2, sqrt2pi, eps_dp, nan_dp, neg_inf_dp
   implicit none
   private

   public :: normal_pdf, normal_cdf, regularized_beta, beta_log_pdf
   public :: safe_log_prob, safe_log1m_prob, clamp01

contains

   elemental real(dp) function clamp01(x) result(y)
      real(dp), intent(in) :: x
      y = min(1.0_dp, max(0.0_dp, x))
   end function clamp01

   elemental real(dp) function normal_pdf(x) result(y)
      real(dp), intent(in) :: x
      y = exp(-0.5_dp*x*x)/sqrt2pi
   end function normal_pdf

   elemental real(dp) function normal_cdf(x) result(y)
      real(dp), intent(in) :: x
      y = 0.5_dp*erfc(-x/sqrt2)
   end function normal_cdf

   elemental real(dp) function safe_log_prob(p) result(v)
      real(dp), intent(in) :: p
      if (p <= 0.0_dp) then
         v = neg_inf_dp()
      else
         v = log(p)
      end if
   end function safe_log_prob

   elemental real(dp) function safe_log1m_prob(p) result(v)
      real(dp), intent(in) :: p
      if (p >= 1.0_dp) then
         v = neg_inf_dp()
      else if (p <= 0.0_dp) then
         v = 0.0_dp
      else
         v = log(1.0_dp-p)
      end if
   end function safe_log1m_prob

   elemental real(dp) function beta_log_pdf(x, a, b) result(v)
      real(dp), intent(in) :: x, a, b
      if (a <= 0.0_dp .or. b <= 0.0_dp) then
         v = nan_dp()
      else if (x <= 0.0_dp .or. x >= 1.0_dp) then
         ! The upstream R density functions explicitly return zero at and
         ! outside the endpoints, even when log=TRUE. Baseline wrappers
         ! reproduce that behavior; this low-level routine uses -infinity.
         v = neg_inf_dp()
      else
         v = (a-1.0_dp)*log(x) + (b-1.0_dp)*log(1.0_dp-x) &
             + log_gamma(a+b) - log_gamma(a) - log_gamma(b)
      end if
   end function beta_log_pdf

   pure real(dp) function beta_cf(a, b, x) result(cf)
      real(dp), intent(in) :: a, b, x
      integer, parameter :: maxit = 100000
      real(dp), parameter :: fpmin = tiny(1.0_dp)/eps_dp
      real(dp), parameter :: tol = 16.0_dp*eps_dp
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

      do m = 1, maxit
         m2 = 2*m
         aa = real(m,dp)*(b-real(m,dp))*x &
              / ((qam+real(m2,dp))*(a+real(m2,dp)))
         d = 1.0_dp+aa*d
         if (abs(d) < fpmin) d = fpmin
         c = 1.0_dp+aa/c
         if (abs(c) < fpmin) c = fpmin
         d = 1.0_dp/d
         h = h*d*c

         aa = -(a+real(m,dp))*(qab+real(m,dp))*x &
              / ((a+real(m2,dp))*(qap+real(m2,dp)))
         d = 1.0_dp+aa*d
         if (abs(d) < fpmin) d = fpmin
         c = 1.0_dp+aa/c
         if (abs(c) < fpmin) c = fpmin
         d = 1.0_dp/d
         del = d*c
         h = h*del
         if (abs(del-1.0_dp) <= tol) exit
      end do
      cf = h
   end function beta_cf

   pure real(dp) function regularized_beta(x, a, b) result(p)
      real(dp), intent(in) :: x, a, b
      real(dp) :: bt

      if (a <= 0.0_dp .or. b <= 0.0_dp) then
         p = nan_dp()
         return
      end if
      if (x <= 0.0_dp) then
         p = 0.0_dp
         return
      else if (x >= 1.0_dp) then
         p = 1.0_dp
         return
      end if

      bt = exp(log_gamma(a+b)-log_gamma(a)-log_gamma(b) &
               + a*log(x)+b*log(1.0_dp-x))
      if (x < (a+1.0_dp)/(a+b+2.0_dp)) then
         p = bt*beta_cf(a,b,x)/a
      else
         p = 1.0_dp-bt*beta_cf(b,a,1.0_dp-x)/b
      end if
      p = clamp01(p)
   end function regularized_beta

end module skewunit_special
