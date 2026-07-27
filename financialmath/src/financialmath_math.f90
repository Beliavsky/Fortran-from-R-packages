! SPDX-License-Identifier: GPL-2.0-only
! Based on FinancialMath 0.1.1, Copyright (C) 2016 Kameron Penn and Jack Schmidt.
module financialmath_math
   use financialmath_kinds, only : dp
   implicit none
   private
   real(dp), parameter :: sqrt_two = sqrt(2.0_dp)
   real(dp), parameter :: inv_sqrt_2pi = 0.3989422804014327_dp
   public :: normal_cdf, normal_pdf, solve_root, finite_real, nearly_equal

   abstract interface
      function scalar_function(x) result(y)
         import dp
         real(dp), intent(in) :: x
         real(dp) :: y
      end function scalar_function
   end interface

contains

   pure real(dp) function normal_cdf(x) result(p)
      real(dp), intent(in) :: x
      p = 0.5_dp*erfc(-x/sqrt_two)
   end function normal_cdf

   pure real(dp) function normal_pdf(x) result(p)
      real(dp), intent(in) :: x
      p = inv_sqrt_2pi*exp(-0.5_dp*x*x)
   end function normal_pdf

   pure logical function finite_real(x) result(ok)
      use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
      real(dp), intent(in) :: x
      ok = ieee_is_finite(x)
   end function finite_real

   pure logical function nearly_equal(a, b, atol, rtol) result(ok)
      real(dp), intent(in) :: a, b
      real(dp), intent(in), optional :: atol, rtol
      real(dp) :: aa, rr
      aa = 1.0e-12_dp
      rr = 1.0e-10_dp
      if (present(atol)) aa = atol
      if (present(rtol)) rr = rtol
      ok = abs(a-b) <= aa + rr*max(abs(a), abs(b))
   end function nearly_equal

   function solve_root(fun, lower, upper, ok, tol, max_iter) result(root)
      procedure(scalar_function) :: fun
      real(dp), intent(in) :: lower, upper
      logical, intent(out) :: ok
      real(dp), intent(in), optional :: tol
      integer, intent(in), optional :: max_iter
      real(dp) :: root
      real(dp) :: a, b, c, fa, fb, fc, eps
      integer :: iter, nmax

      eps = 1.0e-12_dp
      nmax = 200
      if (present(tol)) eps = tol
      if (present(max_iter)) nmax = max_iter
      a = lower
      b = upper
      fa = fun(a)
      fb = fun(b)
      ok = finite_real(fa) .and. finite_real(fb)
      if (.not. ok) then
         root = 0.0_dp
         return
      end if
      if (abs(fa) <= eps) then
         root = a
         return
      end if
      if (abs(fb) <= eps) then
         root = b
         return
      end if
      if (fa*fb > 0.0_dp) then
         ok = .false.
         root = 0.0_dp
         return
      end if
      do iter = 1, nmax
         c = 0.5_dp*(a+b)
         fc = fun(c)
         if (.not. finite_real(fc)) then
            ok = .false.
            root = c
            return
         end if
         if (abs(fc) <= eps .or. abs(b-a) <= eps*(1.0_dp+abs(c))) then
            root = c
            return
         end if
         if (fa*fc <= 0.0_dp) then
            b = c
            fb = fc
         else
            a = c
            fa = fc
         end if
      end do
      root = 0.5_dp*(a+b)
   end function solve_root

end module financialmath_math
