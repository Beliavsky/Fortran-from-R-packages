! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from Risk 1.0 by Saralees Nadarajah and Stephen Chan.
! Copyright (c) 2017 Saralees Nadarajah and Stephen Chan.
module risk_math
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, ieee_is_finite
   use risk_kinds, only : dp
   implicit none
   private

   real(dp), parameter, public :: pi = acos(-1.0_dp)
   real(dp), parameter, public :: sqrt_two = sqrt(2.0_dp)
   real(dp), parameter, public :: sqrt_two_pi = sqrt(2.0_dp*pi)

   public :: quiet_nan, normal_pdf_std, normal_cdf_std, normal_quantile_std
   public :: regularized_beta, student_t_cdf_std

contains

   pure function quiet_nan() result(x)
      real(dp) :: x
      x = ieee_value(0.0_dp, ieee_quiet_nan)
   end function quiet_nan

   pure function normal_pdf_std(x) result(y)
      real(dp), intent(in) :: x
      real(dp) :: y
      y = exp(-0.5_dp*x*x)/sqrt_two_pi
   end function normal_pdf_std

   pure function normal_cdf_std(x) result(p)
      real(dp), intent(in) :: x
      real(dp) :: p
      p = 0.5_dp*erfc(-x/sqrt_two)
   end function normal_cdf_std

   pure function normal_quantile_std(p) result(x)
      real(dp), intent(in) :: p
      real(dp) :: x
      real(dp) :: q, r
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
      real(dp), parameter :: p_low = 0.02425_dp
      real(dp), parameter :: p_high = 1.0_dp - p_low

      if (p < 0.0_dp .or. p > 1.0_dp) then
         x = quiet_nan()
      else if (p <= 0.0_dp) then
         x = -huge(1.0_dp)
      else if (p >= 1.0_dp) then
         x = huge(1.0_dp)
      else if (p < p_low) then
         q = sqrt(-2.0_dp*log(p))
         x = (((((c(1)*q+c(2))*q+c(3))*q+c(4))*q+c(5))*q+c(6))/ &
             ((((d(1)*q+d(2))*q+d(3))*q+d(4))*q+1.0_dp)
      else if (p <= p_high) then
         q = p - 0.5_dp
         r = q*q
         x = (((((a(1)*r+a(2))*r+a(3))*r+a(4))*r+a(5))*r+a(6))*q/ &
             (((((b(1)*r+b(2))*r+b(3))*r+b(4))*r+b(5))*r+1.0_dp)
      else
         q = sqrt(-2.0_dp*log(1.0_dp-p))
         x = -(((((c(1)*q+c(2))*q+c(3))*q+c(4))*q+c(5))*q+c(6))/ &
              ((((d(1)*q+d(2))*q+d(3))*q+d(4))*q+1.0_dp)
      end if

      if (p > 0.0_dp .and. p < 1.0_dp .and. ieee_is_finite(x)) then
         x = x - (normal_cdf_std(x)-p)/normal_pdf_std(x)
      end if
   end function normal_quantile_std

   pure recursive function log_gamma_lanczos(x) result(value)
      real(dp), intent(in) :: x
      real(dp) :: value
      real(dp) :: z, t, s
      integer :: i
      real(dp), parameter :: coeff(9) = [ &
         0.99999999999980993_dp, 676.5203681218851_dp, &
        -1259.1392167224028_dp, 771.32342877765313_dp, &
        -176.61502916214059_dp, 12.507343278686905_dp, &
        -0.13857109526572012_dp, 9.9843695780195716e-6_dp, &
         1.5056327351493116e-7_dp ]

      if (x < 0.5_dp) then
         value = log(pi) - log(sin(pi*x)) - log_gamma_lanczos(1.0_dp-x)
         return
      end if
      z = x - 1.0_dp
      s = coeff(1)
      do i = 2, size(coeff)
         s = s + coeff(i)/(z+real(i-1,dp))
      end do
      t = z + 7.5_dp
      value = 0.5_dp*log(2.0_dp*pi) + (z+0.5_dp)*log(t) - t + log(s)
   end function log_gamma_lanczos

   pure function beta_continued_fraction(a, b, x) result(cf)
      real(dp), intent(in) :: a, b, x
      real(dp) :: cf
      integer, parameter :: max_iter = 300
      real(dp), parameter :: eps = 3.0e-14_dp
      real(dp), parameter :: fpmin = 1.0e-300_dp
      integer :: m, m2
      real(dp) :: aa, c, d, delta, h, qab, qam, qap

      qab = a + b
      qap = a + 1.0_dp
      qam = a - 1.0_dp
      c = 1.0_dp
      d = 1.0_dp - qab*x/qap
      if (abs(d) < fpmin) d = fpmin
      d = 1.0_dp/d
      h = d
      do m = 1, max_iter
         m2 = 2*m
         aa = real(m,dp)*(b-real(m,dp))*x/ &
              ((qam+real(m2,dp))*(a+real(m2,dp)))
         d = 1.0_dp + aa*d
         if (abs(d) < fpmin) d = fpmin
         c = 1.0_dp + aa/c
         if (abs(c) < fpmin) c = fpmin
         d = 1.0_dp/d
         h = h*d*c

         aa = -(a+real(m,dp))*(qab+real(m,dp))*x/ &
              ((a+real(m2,dp))*(qap+real(m2,dp)))
         d = 1.0_dp + aa*d
         if (abs(d) < fpmin) d = fpmin
         c = 1.0_dp + aa/c
         if (abs(c) < fpmin) c = fpmin
         d = 1.0_dp/d
         delta = d*c
         h = h*delta
         if (abs(delta-1.0_dp) <= eps) exit
      end do
      cf = h
   end function beta_continued_fraction

   pure function regularized_beta(x, a, b) result(value)
      real(dp), intent(in) :: x, a, b
      real(dp) :: value
      real(dp) :: front

      if (a <= 0.0_dp .or. b <= 0.0_dp .or. x < 0.0_dp .or. x > 1.0_dp) then
         value = quiet_nan()
      else if (x <= 0.0_dp) then
         value = 0.0_dp
      else if (x >= 1.0_dp) then
         value = 1.0_dp
      else
         front = exp(log_gamma_lanczos(a+b)-log_gamma_lanczos(a)- &
                     log_gamma_lanczos(b)+a*log(x)+b*log(1.0_dp-x))
         if (x < (a+1.0_dp)/(a+b+2.0_dp)) then
            value = front*beta_continued_fraction(a,b,x)/a
         else
            value = 1.0_dp-front*beta_continued_fraction(b,a,1.0_dp-x)/b
         end if
         value = min(1.0_dp,max(0.0_dp,value))
      end if
   end function regularized_beta

   pure function student_t_cdf_std(x, nu) result(p)
      real(dp), intent(in) :: x, nu
      real(dp) :: p, z, ib

      if (nu <= 0.0_dp) then
         p = quiet_nan()
         return
      end if
      if (abs(x) <= tiny(1.0_dp)) then
         p = 0.5_dp
         return
      end if
      z = nu/(nu+x*x)
      ib = regularized_beta(z,0.5_dp*nu,0.5_dp)
      if (x > 0.0_dp) then
         p = 1.0_dp-0.5_dp*ib
      else
         p = 0.5_dp*ib
      end if
   end function student_t_cdf_std

end module risk_math
