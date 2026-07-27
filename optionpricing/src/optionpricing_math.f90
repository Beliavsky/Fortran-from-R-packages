! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
! Derived from OptionPricing 0.1.2 by Wolfgang Hormann and Kemal Dingec.
module optionpricing_math
   use optionpricing_kinds, only : dp, pi
   implicit none
   private
   public :: normal_pdf, normal_cdf, normal_quantile, mean_value, sample_sd
   public :: adaptive_integral, bisection_root, clamp_probability

   abstract interface
      function scalar_function(x) result(y)
         import dp
         real(dp), intent(in) :: x
         real(dp) :: y
      end function scalar_function
   end interface

contains

   pure elemental function normal_pdf(x) result(y)
      real(dp), intent(in) :: x
      real(dp) :: y
      y = exp(-0.5_dp*x*x)/sqrt(2.0_dp*pi)
   end function normal_pdf

   pure elemental function normal_cdf(x) result(y)
      real(dp), intent(in) :: x
      real(dp) :: y
      y = 0.5_dp*erfc(-x/sqrt(2.0_dp))
   end function normal_cdf

   pure elemental function clamp_probability(p) result(q)
      real(dp), intent(in) :: p
      real(dp) :: q
      q = min(max(p, epsilon(1.0_dp)), 1.0_dp-epsilon(1.0_dp))
   end function clamp_probability

   pure elemental function normal_quantile(p) result(x)
      ! Acklam's inverse-normal approximation with one Halley refinement.
      real(dp), intent(in) :: p
      real(dp) :: x
      real(dp), parameter :: a(6) = [ &
         -3.969683028665376e1_dp, 2.209460984245205e2_dp, &
         -2.759285104469687e2_dp, 1.383577518672690e2_dp, &
         -3.066479806614716e1_dp, 2.506628277459239_dp]
      real(dp), parameter :: b(5) = [ &
         -5.447609879822406e1_dp, 1.615858368580409e2_dp, &
         -1.556989798598866e2_dp, 6.680131188771972e1_dp, &
         -1.328068155288572e1_dp]
      real(dp), parameter :: c(6) = [ &
         -7.784894002430293e-3_dp, -3.223964580411365e-1_dp, &
         -2.400758277161838_dp, -2.549732539343734_dp, &
          4.374664141464968_dp, 2.938163982698783_dp]
      real(dp), parameter :: d(4) = [ &
          7.784695709041462e-3_dp, 3.224671290700398e-1_dp, &
          2.445134137142996_dp, 3.754408661907416_dp]
      real(dp), parameter :: plow = 0.02425_dp
      real(dp), parameter :: phigh = 1.0_dp-plow
      real(dp) :: q, r, e, u, pp

      pp = clamp_probability(p)
      if (pp < plow) then
         q = sqrt(-2.0_dp*log(pp))
         x = (((((c(1)*q+c(2))*q+c(3))*q+c(4))*q+c(5))*q+c(6))/ &
             ((((d(1)*q+d(2))*q+d(3))*q+d(4))*q+1.0_dp)
      else if (pp <= phigh) then
         q = pp-0.5_dp
         r = q*q
         x = (((((a(1)*r+a(2))*r+a(3))*r+a(4))*r+a(5))*r+a(6))*q/ &
             (((((b(1)*r+b(2))*r+b(3))*r+b(4))*r+b(5))*r+1.0_dp)
      else
         q = sqrt(-2.0_dp*log(1.0_dp-pp))
         x = -(((((c(1)*q+c(2))*q+c(3))*q+c(4))*q+c(5))*q+c(6))/ &
              ((((d(1)*q+d(2))*q+d(3))*q+d(4))*q+1.0_dp)
      end if
      e = normal_cdf(x)-pp
      u = e/normal_pdf(x)
      x = x-u/(1.0_dp+0.5_dp*x*u)
   end function normal_quantile

   pure function mean_value(x) result(m)
      real(dp), intent(in) :: x(:)
      real(dp) :: m
      if (size(x) == 0) then
         m = 0.0_dp
      else
         m = sum(x)/real(size(x),dp)
      end if
   end function mean_value

   pure function sample_sd(x) result(s)
      real(dp), intent(in) :: x(:)
      real(dp) :: s, m
      if (size(x) <= 1) then
         s = 0.0_dp
      else
         m = mean_value(x)
         s = sqrt(max(0.0_dp, sum((x-m)**2)/real(size(x)-1,dp)))
      end if
   end function sample_sd

   function adaptive_integral(f, a, b, tol, max_depth) result(value)
      procedure(scalar_function) :: f
      real(dp), intent(in) :: a, b
      real(dp), intent(in), optional :: tol
      integer, intent(in), optional :: max_depth
      real(dp) :: value, eps, fa, fb, fm, whole
      integer :: depth

      eps = 1.0e-10_dp
      if (present(tol)) eps = tol
      depth = 24
      if (present(max_depth)) depth = max_depth
      fa = f(a)
      fb = f(b)
      fm = f(0.5_dp*(a+b))
      whole = (b-a)*(fa+4.0_dp*fm+fb)/6.0_dp
      value = recurse(a,b,fa,fm,fb,whole,eps,depth)
   contains
      recursive function recurse(xa, xb, fxa, fxm, fxb, old, local_tol, remaining) result(ans)
         real(dp), intent(in) :: xa, xb, fxa, fxm, fxb, old, local_tol
         integer, intent(in) :: remaining
         real(dp) :: ans, xm, xl, xr, fxl, fxr, left, right, delta
         xm = 0.5_dp*(xa+xb)
         xl = 0.5_dp*(xa+xm)
         xr = 0.5_dp*(xm+xb)
         fxl = f(xl)
         fxr = f(xr)
         left = (xm-xa)*(fxa+4.0_dp*fxl+fxm)/6.0_dp
         right = (xb-xm)*(fxm+4.0_dp*fxr+fxb)/6.0_dp
         delta = left+right-old
         if (remaining <= 0 .or. abs(delta) <= 15.0_dp*local_tol) then
            ans = left+right+delta/15.0_dp
         else
            ans = recurse(xa,xm,fxa,fxl,fxm,left,0.5_dp*local_tol,remaining-1) + &
                  recurse(xm,xb,fxm,fxr,fxb,right,0.5_dp*local_tol,remaining-1)
         end if
      end function recurse
   end function adaptive_integral

   function bisection_root(f, left, right, tol, maxiter, status) result(root)
      procedure(scalar_function) :: f
      real(dp), intent(in) :: left, right
      real(dp), intent(in), optional :: tol
      integer, intent(in), optional :: maxiter
      integer, intent(out), optional :: status
      real(dp) :: root, a, b, fa, fb, fm, eps
      integer :: i, niter, istat
      eps = 1.0e-12_dp
      if (present(tol)) eps = tol
      niter = 200
      if (present(maxiter)) niter = maxiter
      root = 0.5_dp*(left+right)
      a = left
      b = right
      fa = f(a)
      fb = f(b)
      istat = 0
      if (abs(fa) <= eps) then
         root = a
      else if (abs(fb) <= eps) then
         root = b
      else if (fa*fb > 0.0_dp) then
         root = 0.5_dp*(a+b)
         istat = 1
      else
         do i=1,niter
            root = 0.5_dp*(a+b)
            fm = f(root)
            if (abs(fm) <= eps .or. 0.5_dp*(b-a) <= eps) exit
            if (fa*fm <= 0.0_dp) then
               b = root
               fb = fm
            else
               a = root
               fa = fm
            end if
         end do
         if (i > niter) istat = 2
      end if
      if (present(status)) status = istat
   end function bisection_root
end module optionpricing_math
