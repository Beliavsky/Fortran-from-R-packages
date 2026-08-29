! Numerical helpers genuinely missing from the supplied r_compat module.
! SPDX-License-Identifier: GPL-2.0-or-later
module stabledist_numerics
   use r_compat, only: dp, integrate, integrate_result_t
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   implicit none
   private
   public :: bisect_root, golden_max, integrate_split

   abstract interface
      function scalar_fun(x) result(y)
         import dp
         real(dp), intent(in) :: x
         real(dp) :: y
      end function scalar_fun
   end interface

contains

   function bisect_root(fun, a, b, tol, maxiter, ok) result(root)
      procedure(scalar_fun) :: fun
      real(dp), intent(in) :: a, b, tol
      integer, intent(in) :: maxiter
      logical, intent(out), optional :: ok
      real(dp) :: root
      real(dp) :: lo, hi, mid, flo, fhi, fm
      integer :: iter
      logical :: good

      lo = min(a,b)
      hi = max(a,b)
      flo = fun(lo)
      fhi = fun(hi)
      good = ieee_is_finite(flo) .and. ieee_is_finite(fhi)
      if (good) good = (flo == 0.0_dp .or. fhi == 0.0_dp .or. flo*fhi <= 0.0_dp)
      if (.not. good) then
         root = 0.5_dp*(lo+hi)
         if (present(ok)) ok = .false.
         return
      end if
      if (flo == 0.0_dp) then
         root = lo
         if (present(ok)) ok = .true.
         return
      end if
      if (fhi == 0.0_dp) then
         root = hi
         if (present(ok)) ok = .true.
         return
      end if
      do iter = 1, maxiter
         mid = 0.5_dp*(lo+hi)
         fm = fun(mid)
         if (.not. ieee_is_finite(fm)) then
            ! Move away from a singular midpoint without losing the bracket.
            mid = 0.5_dp*(mid+lo)
            fm = fun(mid)
         end if
         if (abs(hi-lo) <= tol*max(1.0_dp,abs(mid)) .or. fm == 0.0_dp) exit
         if (flo*fm <= 0.0_dp) then
            hi = mid
            fhi = fm
         else
            lo = mid
            flo = fm
         end if
      end do
      root = 0.5_dp*(lo+hi)
      if (present(ok)) ok = .true.
   end function bisect_root

   function golden_max(fun, a, b, tol, maxiter) result(xmax)
      procedure(scalar_fun) :: fun
      real(dp), intent(in) :: a, b, tol
      integer, intent(in) :: maxiter
      real(dp) :: xmax
      real(dp), parameter :: gr = 0.6180339887498948482045868343656381_dp
      real(dp) :: lo, hi, x1, x2, f1, f2
      integer :: iter

      lo = min(a,b)
      hi = max(a,b)
      x1 = hi - gr*(hi-lo)
      x2 = lo + gr*(hi-lo)
      f1 = fun(x1)
      f2 = fun(x2)
      do iter = 1, maxiter
         if (abs(hi-lo) <= tol*max(1.0_dp,abs(0.5_dp*(lo+hi)))) exit
         if (f1 < f2) then
            lo = x1
            x1 = x2
            f1 = f2
            x2 = lo + gr*(hi-lo)
            f2 = fun(x2)
         else
            hi = x2
            x2 = x1
            f2 = f1
            x1 = hi - gr*(hi-lo)
            f1 = fun(x1)
         end if
      end do
      xmax = 0.5_dp*(lo+hi)
   end function golden_max

   function integrate_split(fun, a, b, tol, subdivisions, pieces) result(value)
      procedure(scalar_fun) :: fun
      real(dp), intent(in) :: a, b, tol
      integer, intent(in) :: subdivisions, pieces
      real(dp) :: value
      real(dp) :: x0, x1
      integer :: j, np, ns
      type(integrate_result_t) :: ir
      np=max(1,pieces)
      ns=max(8,subdivisions/max(1,np))
      if(mod(ns,2)/=0)ns=ns+1
      value=0.0_dp
      do j=1,np
         x0=a+(b-a)*real(j-1,dp)/real(np,dp)
         x1=a+(b-a)*real(j,dp)/real(np,dp)
         ir=integrate(fun,x0,x1,rel_tol=tol,subdivisions=ns)
         value=value+ir%value
      end do
   end function integrate_split

end module stabledist_numerics
