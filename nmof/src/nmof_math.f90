! SPDX-License-Identifier: GPL-3.0-only
module nmof_math
   use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
   use nmof_kinds, only: dp, pi
   implicit none
   private
   public :: normal_pdf, normal_cdf, normal_quantile
   public :: sort_real, quantile_type7, median_value, standard_deviation
   public :: bisection_root, brent_root, adaptive_simpson
   public :: log_factorial, factorial_real, clamp

   abstract interface
      function scalar_function(x) result(f)
         import dp
         real(dp), intent(in) :: x
         real(dp) :: f
      end function scalar_function
   end interface
contains
   pure elemental real(dp) function clamp(x, lo, hi) result(y)
      real(dp), intent(in) :: x, lo, hi
      y = min(max(x, lo), hi)
   end function clamp

   pure elemental real(dp) function normal_pdf(x) result(f)
      real(dp), intent(in) :: x
      f = exp(-0.5_dp * x * x) / sqrt(2.0_dp * pi)
   end function normal_pdf

   pure elemental real(dp) function normal_cdf(x) result(p)
      real(dp), intent(in) :: x
      p = 0.5_dp * erfc(-x / sqrt(2.0_dp))
   end function normal_cdf

   pure elemental real(dp) function normal_quantile(p) result(x)
      real(dp), intent(in) :: p
      real(dp), parameter :: a1=-3.969683028665376e+01_dp, a2=2.209460984245205e+02_dp
      real(dp), parameter :: a3=-2.759285104469687e+02_dp, a4=1.383577518672690e+02_dp
      real(dp), parameter :: a5=-3.066479806614716e+01_dp, a6=2.506628277459239e+00_dp
      real(dp), parameter :: b1=-5.447609879822406e+01_dp, b2=1.615858368580409e+02_dp
      real(dp), parameter :: b3=-1.556989798598866e+02_dp, b4=6.680131188771972e+01_dp
      real(dp), parameter :: b5=-1.328068155288572e+01_dp
      real(dp), parameter :: c1=-7.784894002430293e-03_dp, c2=-3.223964580411365e-01_dp
      real(dp), parameter :: c3=-2.400758277161838e+00_dp, c4=-2.549732539343734e+00_dp
      real(dp), parameter :: c5=4.374664141464968e+00_dp, c6=2.938163982698783e+00_dp
      real(dp), parameter :: d1=7.784695709041462e-03_dp, d2=3.224671290700398e-01_dp
      real(dp), parameter :: d3=2.445134137142996e+00_dp, d4=3.754408661907416e+00_dp
      real(dp), parameter :: plow=0.02425_dp, phigh=1.0_dp-plow
      real(dp) :: q, r
      if (p <= 0.0_dp) then
         x = -huge(1.0_dp)
      else if (p >= 1.0_dp) then
         x = huge(1.0_dp)
      else if (p < plow) then
         q = sqrt(-2.0_dp*log(p))
         x = (((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6) / ((((d1*q+d2)*q+d3)*q+d4)*q+1.0_dp)
      else if (p <= phigh) then
         q = p-0.5_dp
         r = q*q
         x = (((((a1*r+a2)*r+a3)*r+a4)*r+a5)*r+a6)*q / (((((b1*r+b2)*r+b3)*r+b4)*r+b5)*r+1.0_dp)
      else
         q = sqrt(-2.0_dp*log(1.0_dp-p))
         x = -(((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6) / ((((d1*q+d2)*q+d3)*q+d4)*q+1.0_dp)
      end if
   end function normal_quantile

   subroutine sort_real(x)
      real(dp), intent(inout) :: x(:)
      integer :: i, j
      real(dp) :: key
      do i = 2, size(x)
         key = x(i)
         j = i - 1
         do while (j >= 1)
            if (x(j) <= key) exit
            x(j+1) = x(j)
            j = j - 1
         end do
         x(j+1) = key
      end do
   end subroutine sort_real

   function quantile_type7(x, p) result(q)
      real(dp), intent(in) :: x(:), p
      real(dp) :: q, h, gamma
      real(dp), allocatable :: y(:)
      integer :: n, j
      n = size(x)
      if (n == 0) then
         q = ieee_value(0.0_dp, ieee_quiet_nan)
         return
      end if
      allocate(y(n)); y = x; call sort_real(y)
      if (p <= 0.0_dp) then
         q = y(1)
      else if (p >= 1.0_dp) then
         q = y(n)
      else
         h = 1.0_dp + real(n-1,dp)*p
         j = int(floor(h))
         gamma = h-real(j,dp)
         if (j >= n) then
            q = y(n)
         else
            q = (1.0_dp-gamma)*y(j) + gamma*y(j+1)
         end if
      end if
   end function quantile_type7

   function median_value(x) result(m)
      real(dp), intent(in) :: x(:)
      real(dp) :: m
      m = quantile_type7(x, 0.5_dp)
   end function median_value

   pure function standard_deviation(x) result(sd)
      real(dp), intent(in) :: x(:)
      real(dp) :: sd, mu
      if (size(x) <= 1) then
         sd = 0.0_dp
      else
         mu = sum(x)/real(size(x),dp)
         sd = sqrt(sum((x-mu)**2)/real(size(x)-1,dp))
      end if
   end function standard_deviation

   subroutine bisection_root(fun, lower, upper, root, info, tol, maxiter)
      procedure(scalar_function) :: fun
      real(dp), intent(in) :: lower, upper
      real(dp), intent(out) :: root
      integer, intent(out) :: info
      real(dp), intent(in), optional :: tol
      integer, intent(in), optional :: maxiter
      real(dp) :: a,b,c,fa,fb,fc,eps
      integer :: i, mit
      eps=1.0e-10_dp; if (present(tol)) eps=tol
      mit=200; if (present(maxiter)) mit=maxiter
      a=lower; b=upper; fa=fun(a); fb=fun(b)
      if (abs(fa) <= eps) then; root=a; info=0; return; end if
      if (abs(fb) <= eps) then; root=b; info=0; return; end if
      if (fa*fb>0.0_dp) then; root=0.5_dp*(a+b); info=1; return; end if
      do i=1,mit
         c=0.5_dp*(a+b); fc=fun(c)
         if (abs(fc)<=eps .or. 0.5_dp*abs(b-a)<=eps) then; root=c; info=0; return; end if
         if (fa*fc<=0.0_dp) then; b=c; fb=fc; else; a=c; fa=fc; end if
      end do
      root=0.5_dp*(a+b); info=2
   end subroutine bisection_root

   subroutine brent_root(fun, lower, upper, root, info, tol, maxiter)
      procedure(scalar_function) :: fun
      real(dp), intent(in) :: lower, upper
      real(dp), intent(out) :: root
      integer, intent(out) :: info
      real(dp), intent(in), optional :: tol
      integer, intent(in), optional :: maxiter
      real(dp) :: a,b,fa,fb,x,fx,mid,eps,width,secant
      integer :: iter,mit
      eps=1.0e-10_dp; if(present(tol)) eps=tol
      mit=200; if(present(maxiter)) mit=maxiter
      a=lower; b=upper
      if(a>=b .or. eps<=0.0_dp .or. mit<1) then
         root=0.5_dp*(a+b); info=1; return
      end if
      fa=fun(a); fb=fun(b)
      if(abs(fa)<=eps) then; root=a; info=0; return; end if
      if(abs(fb)<=eps) then; root=b; info=0; return; end if
      if((fa<0.0_dp .and. fb<0.0_dp) .or. (fa>0.0_dp .and. fb>0.0_dp)) then
         root=b; info=1; return
      end if
      do iter=1,mit
         width=b-a; mid=0.5_dp*(a+b)
         if(abs(fb-fa)>tiny(1.0_dp)) then
            secant=b-fb*width/(fb-fa)
         else
            secant=mid
         end if
         if(secant<=a+0.1_dp*width .or. secant>=b-0.1_dp*width) then
            x=mid
         else
            x=secant
         end if
         fx=fun(x)
         if(abs(fx)<=eps .or. 0.5_dp*width<=eps*(1.0_dp+abs(x))) then
            root=x; info=0; return
         end if
         if((fa<0.0_dp .and. fx>0.0_dp) .or. (fa>0.0_dp .and. fx<0.0_dp)) then
            b=x; fb=fx
         else
            a=x; fa=fx
         end if
      end do
      root=0.5_dp*(a+b); info=2
   end subroutine brent_root

   recursive function adaptive_simpson(fun, a, b, tol, max_depth) result(value)
      procedure(scalar_function) :: fun
      real(dp), intent(in) :: a,b,tol
      integer, intent(in), optional :: max_depth
      real(dp) :: value, fa,fb,fc,s
      integer :: depth
      depth=20; if (present(max_depth)) depth=max_depth
      fa=fun(a); fb=fun(b); fc=fun(0.5_dp*(a+b))
      s=(b-a)*(fa+4.0_dp*fc+fb)/6.0_dp
      value=recurse(a,b,fa,fb,fc,s,tol,depth)
   contains
      recursive function recurse(x0,x1,f0,f1,fm,whole,eps,lev) result(ans)
         real(dp), intent(in) :: x0,x1,f0,f1,fm,whole,eps
         integer, intent(in) :: lev
         real(dp) :: ans,xm,xl,xr,fl,fr,left,right
         xm=0.5_dp*(x0+x1); xl=0.5_dp*(x0+xm); xr=0.5_dp*(xm+x1)
         fl=fun(xl); fr=fun(xr)
         left=(xm-x0)*(f0+4.0_dp*fl+fm)/6.0_dp
         right=(x1-xm)*(fm+4.0_dp*fr+f1)/6.0_dp
         if (lev<=0 .or. abs(left+right-whole)<=15.0_dp*eps) then
            ans=left+right+(left+right-whole)/15.0_dp
         else
            ans=recurse(x0,xm,f0,fm,fl,left,0.5_dp*eps,lev-1)+ &
                recurse(xm,x1,fm,f1,fr,right,0.5_dp*eps,lev-1)
         end if
      end function recurse
   end function adaptive_simpson

   pure real(dp) function log_factorial(n) result(v)
      integer, intent(in) :: n
      integer :: i
      v=0.0_dp
      do i=2,n; v=v+log(real(i,dp)); end do
   end function log_factorial

   pure real(dp) function factorial_real(n) result(v)
      integer, intent(in) :: n
      integer :: i
      v=1.0_dp
      do i=2,n; v=v*real(i,dp); end do
   end function factorial_real
end module nmof_math
