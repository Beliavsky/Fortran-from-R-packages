! SPDX-License-Identifier: MIT
module zero_one_special
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
   use zero_one_kinds, only : dp
   implicit none
   private
   public :: normal_cdf, normal_pdf, normal_quantile
   public :: beta_pdf, beta_cdf, beta_quantile
   public :: quiet_nan
contains
   elemental real(dp) function quiet_nan() result(x)
      x = ieee_value(0.0_dp, ieee_quiet_nan)
   end function quiet_nan

   elemental real(dp) function normal_pdf(x) result(v)
      real(dp), intent(in) :: x
      v = exp(-0.5_dp*x*x)/sqrt(2.0_dp*acos(-1.0_dp))
   end function normal_pdf

   elemental real(dp) function normal_cdf(x) result(v)
      real(dp), intent(in) :: x
      v = 0.5_dp*erfc(-x/sqrt(2.0_dp))
   end function normal_cdf

   pure real(dp) function normal_quantile(p) result(x)
      real(dp), intent(in) :: p
      real(dp) :: lo, hi, mid, f, pdf
      integer :: it
      if (p <= 0.0_dp) then
         x = -huge(1.0_dp)
         return
      else if (p >= 1.0_dp) then
         x = huge(1.0_dp)
         return
      end if
      lo = -40.0_dp
      hi = 40.0_dp
      do it = 1, 100
         mid = 0.5_dp*(lo+hi)
         if (normal_cdf(mid) < p) then
            lo = mid
         else
            hi = mid
         end if
      end do
      x = 0.5_dp*(lo+hi)
      do it = 1, 4
         f = normal_cdf(x)-p
         pdf = normal_pdf(x)
         if (pdf <= tiny(1.0_dp)) exit
         mid = x-f/pdf
         if (mid <= lo .or. mid >= hi) exit
         x = mid
      end do
   end function normal_quantile

   elemental real(dp) function beta_pdf(x, a, b) result(v)
      real(dp), intent(in) :: x, a, b
      real(dp) :: lp
      if (a <= 0.0_dp .or. b <= 0.0_dp) then
         v = quiet_nan()
      else if (x < 0.0_dp .or. x > 1.0_dp) then
         v = 0.0_dp
      else if (x <= 0.0_dp) then
         if (a < 1.0_dp) then
            v = huge(1.0_dp)
         else if (abs(a-1.0_dp) <= epsilon(1.0_dp)) then
            v = b
         else
            v = 0.0_dp
         end if
      else if (x >= 1.0_dp) then
         if (b < 1.0_dp) then
            v = huge(1.0_dp)
         else if (abs(b-1.0_dp) <= epsilon(1.0_dp)) then
            v = a
         else
            v = 0.0_dp
         end if
      else
         lp = log_gamma(a+b)-log_gamma(a)-log_gamma(b) + (a-1.0_dp)*log(x) + (b-1.0_dp)*log(1.0_dp-x)
         v = exp(lp)
      end if
   end function beta_pdf

   pure real(dp) function beta_cf(a, b, x) result(h)
      real(dp), intent(in) :: a, b, x
      integer, parameter :: maxit = 300
      real(dp), parameter :: eps = 3.0e-14_dp
      real(dp), parameter :: fpmin = 1.0e-300_dp
      real(dp) :: qab, qap, qam, c, d, aa, del
      integer :: m, m2
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
         aa = real(m,dp)*(b-real(m,dp))*x/((qam+real(m2,dp))*(a+real(m2,dp)))
         d = 1.0_dp+aa*d
         if (abs(d) < fpmin) d = fpmin
         c = 1.0_dp+aa/c
         if (abs(c) < fpmin) c = fpmin
         d = 1.0_dp/d
         h = h*d*c
         aa = -(a+real(m,dp))*(qab+real(m,dp))*x/((a+real(m2,dp))*(qap+real(m2,dp)))
         d = 1.0_dp+aa*d
         if (abs(d) < fpmin) d = fpmin
         c = 1.0_dp+aa/c
         if (abs(c) < fpmin) c = fpmin
         d = 1.0_dp/d
         del = d*c
         h = h*del
         if (abs(del-1.0_dp) <= eps) exit
      end do
   end function beta_cf

   pure real(dp) function beta_cdf(x, a, b) result(v)
      real(dp), intent(in) :: x, a, b
      real(dp) :: bt
      if (a <= 0.0_dp .or. b <= 0.0_dp) then
         v = quiet_nan()
         return
      end if
      if (x <= 0.0_dp) then
         v = 0.0_dp
         return
      else if (x >= 1.0_dp) then
         v = 1.0_dp
         return
      end if
      bt = exp(log_gamma(a+b)-log_gamma(a)-log_gamma(b)+a*log(x)+b*log(1.0_dp-x))
      if (x < (a+1.0_dp)/(a+b+2.0_dp)) then
         v = bt*beta_cf(a,b,x)/a
      else
         v = 1.0_dp-bt*beta_cf(b,a,1.0_dp-x)/b
      end if
      v = min(1.0_dp,max(0.0_dp,v))
   end function beta_cdf

   real(dp) function beta_quantile(p, a, b) result(x)
      real(dp), intent(in) :: p, a, b
      real(dp) :: lo, hi, f, den, cand
      integer :: it
      if (p <= 0.0_dp) then
         x = 0.0_dp
         return
      else if (p >= 1.0_dp) then
         x = 1.0_dp
         return
      end if
      if (a <= 0.0_dp .or. b <= 0.0_dp) then
         x = quiet_nan()
         return
      end if
      lo = 0.0_dp
      hi = 1.0_dp
      x = a/(a+b)
      do it = 1, 120
         f = beta_cdf(x,a,b)-p
         if (f < 0.0_dp) then
            lo = x
         else
            hi = x
         end if
         den = beta_pdf(x,a,b)
         if (den > 1.0e-300_dp) then
            cand = x-f/den
         else
            cand = 0.5_dp*(lo+hi)
         end if
         if (cand <= lo .or. cand >= hi) cand = 0.5_dp*(lo+hi)
         x = cand
         if (hi-lo <= 2.0e-14_dp*max(1.0_dp,abs(x))) exit
      end do
   end function beta_quantile
end module zero_one_special
