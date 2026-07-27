! Part of the experimental modern Fortran translation of tseries 0.10-62.
! This file was created or modified for the Fortran project on 2026-07-23.
! Original tseries authors retain copyright; see NOTICE and ORIGIN.md.
! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only

module tseries_special
   use tseries_kinds, only : dp
   implicit none
   private

   public :: normal_cdf
   public :: chi_square_cdf
   public :: f_cdf
   public :: regularized_gamma_p
   public :: regularized_beta

contains

   pure real(dp) function normal_cdf(x) result(p)
      real(dp), intent(in) :: x
      p = 0.5_dp * erfc(-x / sqrt(2.0_dp))
   end function normal_cdf

   recursive pure real(dp) function log_gamma_lanczos(x) result(value)
      real(dp), intent(in) :: x
      real(dp), parameter :: coeff(9) = [ &
         0.99999999999980993_dp, 676.5203681218851_dp, &
        -1259.1392167224028_dp, 771.32342877765313_dp, &
        -176.61502916214059_dp, 12.507343278686905_dp, &
        -0.13857109526572012_dp, 9.9843695780195716e-6_dp, &
         1.5056327351493116e-7_dp ]
      real(dp), parameter :: pi = acos(-1.0_dp)
      real(dp) :: y, t, sumc
      integer :: i

      if (x < 0.5_dp) then
         value = log(pi) - log(sin(pi*x)) - log_gamma_lanczos(1.0_dp-x)
         return
      end if
      y = x - 1.0_dp
      sumc = coeff(1)
      do i = 2, size(coeff)
         sumc = sumc + coeff(i)/(y+real(i-1,dp))
      end do
      t = y + 7.5_dp
      value = 0.5_dp*log(2.0_dp*pi) + (y+0.5_dp)*log(t) - t + log(sumc)
   end function log_gamma_lanczos

   pure real(dp) function regularized_gamma_p(a, x) result(p)
      real(dp), intent(in) :: a, x
      integer, parameter :: max_iter = 10000
      real(dp), parameter :: eps = 2.0e-15_dp
      real(dp), parameter :: tiny = 1.0e-300_dp
      real(dp) :: ap, del, sumv, b, c, d, h, an, q
      integer :: n

      if (a <= 0.0_dp .or. x < 0.0_dp) then
         p = 0.0_dp
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
         do n = 1, max_iter
            ap = ap + 1.0_dp
            del = del*x/ap
            sumv = sumv + del
            if (abs(del) <= abs(sumv)*eps) exit
         end do
         p = sumv*exp(-x+a*log(x)-log_gamma_lanczos(a))
      else
         b = x + 1.0_dp - a
         c = 1.0_dp/tiny
         d = 1.0_dp/b
         h = d
         do n = 1, max_iter
            an = -real(n,dp)*(real(n,dp)-a)
            b = b + 2.0_dp
            d = an*d + b
            if (abs(d) < tiny) d = tiny
            c = b + an/c
            if (abs(c) < tiny) c = tiny
            d = 1.0_dp/d
            del = d*c
            h = h*del
            if (abs(del-1.0_dp) <= eps) exit
         end do
         q = exp(-x+a*log(x)-log_gamma_lanczos(a))*h
         p = 1.0_dp - q
      end if
      p = max(0.0_dp, min(1.0_dp, p))
   end function regularized_gamma_p

   pure real(dp) function beta_continued_fraction(a, b, x) result(cf)
      real(dp), intent(in) :: a, b, x
      integer, parameter :: max_iter = 10000
      real(dp), parameter :: eps = 3.0e-14_dp
      real(dp), parameter :: tiny = 1.0e-300_dp
      real(dp) :: qab, qap, qam, c, d, h, aa, del
      integer :: m, m2

      qab = a+b
      qap = a+1.0_dp
      qam = a-1.0_dp
      c = 1.0_dp
      d = 1.0_dp-qab*x/qap
      if (abs(d) < tiny) d = tiny
      d = 1.0_dp/d
      h = d
      do m = 1, max_iter
         m2 = 2*m
         aa = real(m,dp)*(b-real(m,dp))*x/((qam+real(m2,dp))*(a+real(m2,dp)))
         d = 1.0_dp+aa*d
         if (abs(d) < tiny) d = tiny
         c = 1.0_dp+aa/c
         if (abs(c) < tiny) c = tiny
         d = 1.0_dp/d
         h = h*d*c
         aa = -(a+real(m,dp))*(qab+real(m,dp))*x/((a+real(m2,dp))*(qap+real(m2,dp)))
         d = 1.0_dp+aa*d
         if (abs(d) < tiny) d = tiny
         c = 1.0_dp+aa/c
         if (abs(c) < tiny) c = tiny
         d = 1.0_dp/d
         del = d*c
         h = h*del
         if (abs(del-1.0_dp) <= eps) exit
      end do
      cf = h
   end function beta_continued_fraction

   pure real(dp) function regularized_beta(x, a, b) result(value)
      real(dp), intent(in) :: x, a, b
      real(dp) :: bt

      if (x <= 0.0_dp) then
         value = 0.0_dp
         return
      else if (x >= 1.0_dp) then
         value = 1.0_dp
         return
      else if (a <= 0.0_dp .or. b <= 0.0_dp) then
         value = 0.0_dp
         return
      end if
      bt = exp(log_gamma_lanczos(a+b)-log_gamma_lanczos(a)-log_gamma_lanczos(b) &
           + a*log(x)+b*log(1.0_dp-x))
      if (x < (a+1.0_dp)/(a+b+2.0_dp)) then
         value = bt*beta_continued_fraction(a,b,x)/a
      else
         value = 1.0_dp-bt*beta_continued_fraction(b,a,1.0_dp-x)/b
      end if
      value = max(0.0_dp, min(1.0_dp, value))
   end function regularized_beta

   pure real(dp) function chi_square_cdf(x, df) result(p)
      real(dp), intent(in) :: x, df
      if (x <= 0.0_dp) then
         p = 0.0_dp
      else
         p = regularized_gamma_p(0.5_dp*df, 0.5_dp*x)
      end if
   end function chi_square_cdf

   pure real(dp) function f_cdf(x, df1, df2) result(p)
      real(dp), intent(in) :: x, df1, df2
      real(dp) :: z
      if (x <= 0.0_dp) then
         p = 0.0_dp
      else
         z = (df1*x)/(df1*x+df2)
         p = regularized_beta(z, 0.5_dp*df1, 0.5_dp*df2)
      end if
   end function f_cdf

end module tseries_special
