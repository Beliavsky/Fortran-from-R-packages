! Part of the experimental modern Fortran translation of rugarch 1.5-6.
! This file was created or modified for the Fortran project on 2026-07-23.
! Original rugarch authors retain copyright; see NOTICE and ORIGIN.md.
! SPDX-License-Identifier: GPL-3.0-only

module rugarch_math
   use rugarch_kinds, only : dp
   implicit none
   private

   real(dp), parameter, public :: pi = acos(-1.0_dp)
   real(dp), parameter :: sqrt_two = sqrt(2.0_dp)
   real(dp), parameter :: tiny_prob = 1.0e-14_dp

   public :: normal_pdf, normal_cdf, normal_quantile
   public :: regularized_gamma_p, inverse_regularized_gamma_p
   public :: regularized_beta, student_t_pdf, student_t_cdf
   public :: student_t_quantile, beta_fn

contains

   pure elemental function normal_pdf(x) result(value)
      real(dp), intent(in) :: x
      real(dp) :: value

      value = exp(-0.5_dp*x*x) / sqrt(2.0_dp*pi)
   end function normal_pdf

   pure elemental function normal_cdf(x) result(value)
      real(dp), intent(in) :: x
      real(dp) :: value

      value = 0.5_dp*erfc(-x/sqrt_two)
   end function normal_cdf

   pure elemental function normal_quantile(p) result(x)
      real(dp), intent(in) :: p
      real(dp) :: x, q, r
      real(dp), parameter :: a(6) = [ &
         -3.969683028665376e+01_dp,  2.209460984245205e+02_dp, &
         -2.759285104469687e+02_dp,  1.383577518672690e+02_dp, &
         -3.066479806614716e+01_dp,  2.506628277459239e+00_dp ]
      real(dp), parameter :: b(5) = [ &
         -5.447609879822406e+01_dp,  1.615858368580409e+02_dp, &
         -1.556989798598866e+02_dp,  6.680131188771972e+01_dp, &
         -1.328068155288572e+01_dp ]
      real(dp), parameter :: c(6) = [ &
         -7.784894002430293e-03_dp, -3.223964580411365e-01_dp, &
         -2.400758277161838e+00_dp, -2.549732539343734e+00_dp, &
          4.374664141464968e+00_dp,  2.938163982698783e+00_dp ]
      real(dp), parameter :: d(4) = [ &
          7.784695709041462e-03_dp,  3.224671290700398e-01_dp, &
          2.445134137142996e+00_dp,  3.754408661907416e+00_dp ]
      real(dp), parameter :: plow = 0.02425_dp
      real(dp), parameter :: phigh = 1.0_dp - plow

      if (p <= 0.0_dp) then
         x = -huge(1.0_dp)
      else if (p >= 1.0_dp) then
         x = huge(1.0_dp)
      else if (p < plow) then
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
   end function normal_quantile

   pure elemental function beta_fn(a, b) result(value)
      real(dp), intent(in) :: a, b
      real(dp) :: value

      value = exp(log_gamma(a) + log_gamma(b) - log_gamma(a+b))
   end function beta_fn

   pure function regularized_gamma_p(a, x) result(value)
      real(dp), intent(in) :: a, x
      real(dp) :: value
      real(dp) :: ap, del, sumv, b, c, d, h, an
      integer :: i
      integer, parameter :: max_iter = 500
      real(dp), parameter :: eps = 5.0e-15_dp
      real(dp), parameter :: fpmin = tiny(1.0_dp)/eps

      if (a <= 0.0_dp .or. x < 0.0_dp) then
         value = 0.0_dp
         return
      end if
      if (abs(x) <= tiny(1.0_dp)) then
         value = 0.0_dp
         return
      end if

      if (x < a + 1.0_dp) then
         ap = a
         sumv = 1.0_dp/a
         del = sumv
         do i = 1, max_iter
            ap = ap + 1.0_dp
            del = del*x/ap
            sumv = sumv + del
            if (abs(del) <= abs(sumv)*eps) exit
         end do
         value = sumv*exp(-x + a*log(x) - log_gamma(a))
      else
         b = x + 1.0_dp - a
         c = 1.0_dp/fpmin
         d = 1.0_dp/b
         h = d
         do i = 1, max_iter
            an = -real(i,dp)*(real(i,dp)-a)
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
         value = 1.0_dp - exp(-x + a*log(x) - log_gamma(a))*h
      end if
      value = max(0.0_dp, min(1.0_dp, value))
   end function regularized_gamma_p

   pure function inverse_regularized_gamma_p(a, p) result(x)
      real(dp), intent(in) :: a, p
      real(dp) :: x, lo, hi, candidate, value, density, step, slope, z, base
      integer :: i

      if (p <= 0.0_dp) then
         x = 0.0_dp
         return
      else if (p >= 1.0_dp) then
         x = huge(1.0_dp)
         return
      end if

      ! Start from a lower-tail power approximation or the Wilson-Hilferty
      ! approximation, then use safeguarded Halley steps.  The former fixed
      ! 120-step bisection spent almost all qged time evaluating the gamma CDF.
      if (p < 0.5_dp) then
         x = (p*gamma(a+1.0_dp))**(1.0_dp/a)
      else
         z = normal_quantile(p)
         base = 1.0_dp - 1.0_dp/(9.0_dp*a) + z/(3.0_dp*sqrt(a))
         if (base > 0.0_dp) then
            x = a*base**3
         else
            x = max(0.1_dp, a)
         end if
      end if
      lo = 0.0_dp
      hi = max(1.0_dp, a, 1.25_dp*x)
      do while (regularized_gamma_p(a,hi) < p .and. hi < 1.0e12_dp)
         hi = 2.0_dp*hi
      end do
      x = max(tiny(1.0_dp), min(hi, x))
      do i = 1, 20
         value = regularized_gamma_p(a,x)
         if (value < p) then
            lo = x
         else
            hi = x
         end if
         if (abs(value-p) <= 8.0_dp*epsilon(p)*max(p,1.0_dp-p)) exit
         density = exp((a-1.0_dp)*log(x)-x-log_gamma(a))
         if (density > tiny(1.0_dp)) then
            step = (value-p)/density
            slope = (a-1.0_dp)/x-1.0_dp
            step = step/(1.0_dp-0.5_dp*step*slope)
            candidate = x-step
         else
            candidate = 0.5_dp*(lo+hi)
         end if
         if (.not. (candidate > lo .and. candidate < hi)) candidate=0.5_dp*(lo+hi)
         x = candidate
      end do
      x = max(0.0_dp,x)
   end function inverse_regularized_gamma_p

   pure function beta_cont_frac(a, b, x) result(value)
      real(dp), intent(in) :: a, b, x
      real(dp) :: value
      real(dp) :: qab, qap, qam, c, d, h, aa, del
      integer :: m, m2
      integer, parameter :: max_iter = 500
      real(dp), parameter :: eps = 5.0e-15_dp
      real(dp), parameter :: fpmin = tiny(1.0_dp)/eps

      qab = a+b
      qap = a+1.0_dp
      qam = a-1.0_dp
      c = 1.0_dp
      d = 1.0_dp - qab*x/qap
      if (abs(d) < fpmin) d = fpmin
      d = 1.0_dp/d
      h = d
      do m = 1, max_iter
         m2 = 2*m
         aa = real(m,dp)*(b-real(m,dp))*x / &
              ((qam+real(m2,dp))*(a+real(m2,dp)))
         d = 1.0_dp + aa*d
         if (abs(d) < fpmin) d = fpmin
         c = 1.0_dp + aa/c
         if (abs(c) < fpmin) c = fpmin
         d = 1.0_dp/d
         h = h*d*c
         aa = -(a+real(m,dp))*(qab+real(m,dp))*x / &
              ((a+real(m2,dp))*(qap+real(m2,dp)))
         d = 1.0_dp + aa*d
         if (abs(d) < fpmin) d = fpmin
         c = 1.0_dp + aa/c
         if (abs(c) < fpmin) c = fpmin
         d = 1.0_dp/d
         del = d*c
         h = h*del
         if (abs(del-1.0_dp) <= eps) exit
      end do
      value = h
   end function beta_cont_frac

   pure function regularized_beta(x, a, b) result(value)
      real(dp), intent(in) :: x, a, b
      real(dp) :: value, bt

      if (x <= 0.0_dp) then
         value = 0.0_dp
         return
      else if (x >= 1.0_dp) then
         value = 1.0_dp
         return
      end if
      bt = exp(log_gamma(a+b)-log_gamma(a)-log_gamma(b) + &
               a*log(x) + b*log(1.0_dp-x))
      if (x < (a+1.0_dp)/(a+b+2.0_dp)) then
         value = bt*beta_cont_frac(a,b,x)/a
      else
         value = 1.0_dp - bt*beta_cont_frac(b,a,1.0_dp-x)/b
      end if
      value = max(0.0_dp, min(1.0_dp, value))
   end function regularized_beta

   pure elemental function student_t_pdf(x, nu) result(value)
      real(dp), intent(in) :: x, nu
      real(dp) :: value

      if (nu <= 0.0_dp) then
         value = 0.0_dp
      else
         value = exp(log_gamma(0.5_dp*(nu+1.0_dp))-log_gamma(0.5_dp*nu)) / &
                 sqrt(pi*nu) * (1.0_dp+x*x/nu)**(-0.5_dp*(nu+1.0_dp))
      end if
   end function student_t_pdf

   pure elemental function student_t_cdf(x, nu) result(value)
      real(dp), intent(in) :: x, nu
      real(dp) :: value, ib

      if (nu <= 0.0_dp) then
         value = 0.0_dp
         return
      end if
      if (abs(x) <= tiny(1.0_dp)) then
         value = 0.5_dp
      else
         ib = regularized_beta(nu/(nu+x*x), 0.5_dp*nu, 0.5_dp)
         if (x > 0.0_dp) then
            value = 1.0_dp - 0.5_dp*ib
         else
            value = 0.5_dp*ib
         end if
      end if
   end function student_t_cdf

   pure elemental function student_t_quantile(p, nu) result(x)
      real(dp), intent(in) :: p, nu
      real(dp) :: x, lo, hi, candidate, value, density, step, z, z2
      integer :: i

      if (p <= 0.0_dp) then
         x = -huge(1.0_dp)
         return
      else if (p >= 1.0_dp) then
         x = huge(1.0_dp)
         return
      else if (abs(p-0.5_dp) <= epsilon(p)) then
         x = 0.0_dp
         return
      end if

      if (nu <= 0.0_dp) then
         x = huge(1.0_dp)
         return
      end if

      ! A Cornish-Fisher expansion gives a close initial value for ordinary
      ! degrees of freedom.  Safeguarded Newton steps then retain the robust
      ! bracketing of the old 100-evaluation bisection.
      z = normal_quantile(p)
      z2 = z*z
      x = z + z*(z2+1.0_dp)/(4.0_dp*nu) + &
          z*(5.0_dp*z2*z2+16.0_dp*z2+3.0_dp)/(96.0_dp*nu*nu) + &
          z*(3.0_dp*z2**3+19.0_dp*z2*z2+17.0_dp*z2-15.0_dp)/(384.0_dp*nu**3)
      lo = min(-1.0_dp,1.25_dp*x)
      hi = max( 1.0_dp,1.25_dp*x)
      do while (student_t_cdf(lo,nu) > p)
         lo=2.0_dp*lo
      end do
      do while (student_t_cdf(hi,nu) < p)
         hi=2.0_dp*hi
      end do
      x=max(lo,min(hi,x))
      do i = 1, 16
         value=student_t_cdf(x,nu)
         if (value < p) then
            lo=x
         else
            hi=x
         end if
         if (abs(value-p) <= 8.0_dp*epsilon(p)*max(p,1.0_dp-p)) exit
         density=student_t_pdf(x,nu)
         if (density > tiny(1.0_dp)) then
            step=(value-p)/density
            candidate=x-step
         else
            candidate=0.5_dp*(lo+hi)
         end if
         if (.not. (candidate > lo .and. candidate < hi)) candidate=0.5_dp*(lo+hi)
         x=candidate
      end do
   end function student_t_quantile

end module rugarch_math
