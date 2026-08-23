module aep_special
   implicit none
   private
   integer, parameter, public :: dp = kind(1.0d0)
   real(dp), parameter :: eps = epsilon(1.0_dp)
   public :: reg_gamma_p, gamma_quantile, digamma_aep, golden_min, bisect_root
contains
   pure real(dp) function reg_gamma_p(a, x) result(p)
      real(dp), intent(in) :: a, x
      integer, parameter :: maxit = 10000
      real(dp), parameter :: fpmin = tiny(1.0_dp)/eps
      integer :: n
      real(dp) :: ap, del, sumv, b, c, d, h, an, q
      if (a <= 0.0_dp .or. x < 0.0_dp) then
         p = ieee_nan()
         return
      end if
      if (x == 0.0_dp) then
         p = 0.0_dp
         return
      end if
      if (x < a + 1.0_dp) then
         ap = a
         sumv = 1.0_dp/a
         del = sumv
         do n = 1, maxit
            ap = ap + 1.0_dp
            del = del*x/ap
            sumv = sumv + del
            if (abs(del) <= abs(sumv)*8.0_dp*eps) exit
         end do
         p = sumv*exp(-x + a*log(x) - log_gamma(a))
      else
         b = x + 1.0_dp - a
         c = 1.0_dp/fpmin
         d = 1.0_dp/b
         h = d
         do n = 1, maxit
            an = -real(n,dp)*(real(n,dp)-a)
            b = b + 2.0_dp
            d = an*d + b
            if (abs(d) < fpmin) d = fpmin
            c = b + an/c
            if (abs(c) < fpmin) c = fpmin
            d = 1.0_dp/d
            del = d*c
            h = h*del
            if (abs(del - 1.0_dp) <= 8.0_dp*eps) exit
         end do
         q = exp(-x + a*log(x) - log_gamma(a))*h
         p = 1.0_dp - q
      end if
      p = max(0.0_dp, min(1.0_dp, p))
   end function reg_gamma_p

   pure real(dp) function gamma_quantile(prob, a) result(x)
      real(dp), intent(in) :: prob, a
      real(dp) :: lo, hi, mid, pmid
      integer :: it
      if (a <= 0.0_dp .or. prob < 0.0_dp .or. prob > 1.0_dp) then
         x = ieee_nan(); return
      end if
      if (prob == 0.0_dp) then
         x = 0.0_dp; return
      else if (prob == 1.0_dp) then
         x = huge(1.0_dp); return
      end if
      lo = 0.0_dp
      hi = max(1.0_dp, a)
      do while (reg_gamma_p(a, hi) < prob .and. hi < huge(1.0_dp)/4.0_dp)
         hi = 2.0_dp*hi
      end do
      do it = 1, 120
         mid = 0.5_dp*(lo + hi)
         pmid = reg_gamma_p(a, mid)
         if (pmid < prob) then
            lo = mid
         else
            hi = mid
         end if
         if (hi-lo <= 8.0_dp*epsilon(mid)*max(1.0_dp,mid)) exit
      end do
      x = 0.5_dp*(lo+hi)
   end function gamma_quantile

   pure real(dp) function digamma_aep(x) result(y)
      real(dp), intent(in) :: x
      real(dp) :: z, inv, inv2
      if (x <= 0.0_dp) then
         y = ieee_nan(); return
      end if
      z = x; y = 0.0_dp
      do while (z < 8.0_dp)
         y = y - 1.0_dp/z
         z = z + 1.0_dp
      end do
      inv = 1.0_dp/z; inv2 = inv*inv
      y = y + log(z) - 0.5_dp*inv - inv2*(1.0_dp/12.0_dp - inv2*(1.0_dp/120.0_dp - inv2*(1.0_dp/252.0_dp)))
   end function digamma_aep

   real(dp) function golden_min(f, a, b, tol) result(xmin)
      interface
         function f(x) result(y)
            import dp
            real(dp), intent(in) :: x
            real(dp) :: y
         end function f
      end interface
      real(dp), intent(in) :: a, b, tol
      real(dp), parameter :: gr = 0.6180339887498948482_dp
      real(dp) :: lo, hi, x1, x2, f1, f2
      integer :: it
      lo=a; hi=b; x1=hi-gr*(hi-lo); x2=lo+gr*(hi-lo); f1=f(x1); f2=f(x2)
      do it=1,500
         if (abs(hi-lo) <= tol*max(1.0_dp,abs(lo)+abs(hi))) exit
         if (f1 > f2) then
            lo=x1; x1=x2; f1=f2; x2=lo+gr*(hi-lo); f2=f(x2)
         else
            hi=x2; x2=x1; f2=f1; x1=hi-gr*(hi-lo); f1=f(x1)
         end if
      end do
      xmin=0.5_dp*(lo+hi)
   end function golden_min

   real(dp) function bisect_root(f, a, b, tol, ok) result(root)
      interface
         function f(x) result(y)
            import dp
            real(dp), intent(in) :: x
            real(dp) :: y
         end function f
      end interface
      real(dp), intent(in) :: a,b,tol
      logical, intent(out) :: ok
      real(dp) :: lo,hi,mid,fl,fm
      integer :: it
      lo=a; hi=b; fl=f(lo)
      if (fl*f(hi)>0.0_dp) then; ok=.false.; root=0.5_dp*(a+b); return; end if
      do it=1,200
         mid=0.5_dp*(lo+hi); fm=f(mid)
         if (abs(fm)<=tol .or. abs(hi-lo)<=tol) exit
         if (fl*fm<=0.0_dp) then; hi=mid; else; lo=mid; fl=fm; end if
      end do
      root=0.5_dp*(lo+hi); ok=.true.
   end function bisect_root

   pure real(dp) function ieee_nan() result(x)
      use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
      x = ieee_value(0.0_dp, ieee_quiet_nan)
   end function ieee_nan
end module aep_special
