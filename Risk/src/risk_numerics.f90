! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from Risk 1.0 by Saralees Nadarajah and Stephen Chan.
! Copyright (c) 2017 Saralees Nadarajah and Stephen Chan.
module risk_numerics
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use risk_kinds, only : dp
   use risk_math, only : pi, quiet_nan
   implicit none
   private

   integer, parameter :: transform_finite = 0
   integer, parameter :: transform_whole_line = 1
   integer, parameter :: transform_left_infinite = 2
   integer, parameter :: transform_right_infinite = 3

   type, abstract, public :: numeric_function
   contains
      procedure(evaluate_method), deferred :: evaluate
   end type numeric_function

   abstract interface
      function evaluate_method(self, x) result(y)
         import :: numeric_function, dp
         class(numeric_function), intent(in) :: self
         real(dp), intent(in) :: x
         real(dp) :: y
      end function evaluate_method
   end interface

   public :: integrate, solve_bisection
   public :: is_negative_infinity, is_positive_infinity

contains

   pure function is_positive_infinity(x) result(is_inf)
      real(dp), intent(in) :: x
      logical :: is_inf
      is_inf = x >= 0.5_dp*huge(1.0_dp)
   end function is_positive_infinity

   pure function is_negative_infinity(x) result(is_inf)
      real(dp), intent(in) :: x
      logical :: is_inf
      is_inf = x <= -0.5_dp*huge(1.0_dp)
   end function is_negative_infinity

   recursive function integrate(fun, a, b, abs_tol, rel_tol, max_depth) result(value)
      class(numeric_function), intent(in) :: fun
      real(dp), intent(in) :: a, b
      real(dp), intent(in), optional :: abs_tol, rel_tol
      integer, intent(in), optional :: max_depth
      real(dp) :: value
      real(dp) :: atol, rtol
      integer :: depth, transform

      atol = 1.0e-9_dp
      rtol = 1.0e-9_dp
      depth = 30
      if (present(abs_tol)) atol = max(abs_tol,0.0_dp)
      if (present(rel_tol)) rtol = max(rel_tol,0.0_dp)
      if (present(max_depth)) depth = max(max_depth,1)

      if (.not. (a < b) .and. .not. (a > b)) then
         value = 0.0_dp
         return
      else if (a > b) then
         value = -integrate(fun,b,a,atol,rtol,depth)
         return
      end if

      if (is_negative_infinity(a) .and. is_positive_infinity(b)) then
         transform = transform_whole_line
         value = integrate_finite(fun,0.0_dp,1.0_dp,a,b,transform,atol,rtol,depth)
      else if (is_negative_infinity(a)) then
         transform = transform_left_infinite
         value = integrate_finite(fun,0.0_dp,1.0_dp,a,b,transform,atol,rtol,depth)
      else if (is_positive_infinity(b)) then
         transform = transform_right_infinite
         value = integrate_finite(fun,0.0_dp,1.0_dp,a,b,transform,atol,rtol,depth)
      else
         transform = transform_finite
         value = integrate_finite(fun,a,b,a,b,transform,atol,rtol,depth)
      end if
   end function integrate

   function integrate_finite(fun, a, b, original_a, original_b, transform, &
                             abs_tol, rel_tol, max_depth) result(value)
      class(numeric_function), intent(in) :: fun
      real(dp), intent(in) :: a, b, original_a, original_b, abs_tol, rel_tol
      integer, intent(in) :: transform, max_depth
      real(dp) :: value
      real(dp) :: fa, fm, fb, midpoint, whole, tolerance

      midpoint = 0.5_dp*(a+b)
      fa = transformed_value(fun,a,original_a,original_b,transform)
      fm = transformed_value(fun,midpoint,original_a,original_b,transform)
      fb = transformed_value(fun,b,original_a,original_b,transform)
      if (.not. ieee_is_finite(fa) .or. .not. ieee_is_finite(fm) .or. &
          .not. ieee_is_finite(fb)) then
         value = quiet_nan()
         return
      end if
      whole = (b-a)*(fa+4.0_dp*fm+fb)/6.0_dp
      tolerance = max(abs_tol,rel_tol*abs(whole))
      value = adaptive_simpson(fun,a,b,original_a,original_b,transform, &
                               fa,fm,fb,whole,tolerance,max_depth)
   end function integrate_finite

   recursive function adaptive_simpson(fun, a, b, original_a, original_b, &
                                       transform, fa, fm, fb, whole, tol, depth) result(value)
      class(numeric_function), intent(in) :: fun
      real(dp), intent(in) :: a, b, original_a, original_b
      integer, intent(in) :: transform
      real(dp), intent(in) :: fa, fm, fb, whole, tol
      integer, intent(in) :: depth
      real(dp) :: value
      real(dp) :: left_mid, right_mid, center, f_left_mid, f_right_mid
      real(dp) :: left, right, delta

      center = 0.5_dp*(a+b)
      left_mid = 0.5_dp*(a+center)
      right_mid = 0.5_dp*(center+b)
      f_left_mid = transformed_value(fun,left_mid,original_a,original_b,transform)
      f_right_mid = transformed_value(fun,right_mid,original_a,original_b,transform)
      if (.not. ieee_is_finite(f_left_mid) .or. .not. ieee_is_finite(f_right_mid)) then
         value = quiet_nan()
         return
      end if

      left = (center-a)*(fa+4.0_dp*f_left_mid+fm)/6.0_dp
      right = (b-center)*(fm+4.0_dp*f_right_mid+fb)/6.0_dp
      delta = left+right-whole
      if (depth <= 0 .or. abs(delta) <= 15.0_dp*tol) then
         value = left+right+delta/15.0_dp
      else
         value = adaptive_simpson(fun,a,center,original_a,original_b,transform, &
                                  fa,f_left_mid,fm,left,0.5_dp*tol,depth-1)+ &
                 adaptive_simpson(fun,center,b,original_a,original_b,transform, &
                                  fm,f_right_mid,fb,right,0.5_dp*tol,depth-1)
      end if
   end function adaptive_simpson

   function transformed_value(fun, t, a, b, transform) result(y)
      class(numeric_function), intent(in) :: fun
      real(dp), intent(in) :: t, a, b
      integer, intent(in) :: transform
      real(dp) :: y
      real(dp) :: angle, c, one_minus_t, x

      select case (transform)
      case (transform_finite)
         y = fun%evaluate(t)
      case (transform_whole_line)
         if (t <= 0.0_dp .or. t >= 1.0_dp) then
            y = 0.0_dp
         else
            angle = pi*(t-0.5_dp)
            c = cos(angle)
            x = tan(angle)
            y = fun%evaluate(x)*pi/(c*c)
         end if
      case (transform_left_infinite)
         if (t <= 0.0_dp) then
            y = fun%evaluate(b)
         else if (t >= 1.0_dp) then
            y = 0.0_dp
         else
            one_minus_t = 1.0_dp-t
            x = b-t/one_minus_t
            y = fun%evaluate(x)/(one_minus_t*one_minus_t)
         end if
      case (transform_right_infinite)
         if (t <= 0.0_dp) then
            y = fun%evaluate(a)
         else if (t >= 1.0_dp) then
            y = 0.0_dp
         else
            one_minus_t = 1.0_dp-t
            x = a+t/one_minus_t
            y = fun%evaluate(x)/(one_minus_t*one_minus_t)
         end if
      case default
         y = quiet_nan()
      end select
   end function transformed_value

   function solve_bisection(fun, lower, upper, abs_tol, rel_tol, max_iter) result(root)
      class(numeric_function), intent(in) :: fun
      real(dp), intent(in) :: lower, upper
      real(dp), intent(in), optional :: abs_tol, rel_tol
      integer, intent(in), optional :: max_iter
      real(dp) :: root
      real(dp) :: lo, hi, mid, flo, fhi, fmid, atol, rtol
      integer :: i, iterations

      atol = 1.0e-10_dp
      rtol = 1.0e-10_dp
      iterations = 200
      if (present(abs_tol)) atol = max(abs_tol,0.0_dp)
      if (present(rel_tol)) rtol = max(rel_tol,0.0_dp)
      if (present(max_iter)) iterations = max(max_iter,1)

      lo = lower
      hi = upper
      flo = fun%evaluate(lo)
      fhi = fun%evaluate(hi)
      if (.not. ieee_is_finite(flo) .or. .not. ieee_is_finite(fhi) .or. flo*fhi > 0.0_dp) then
         root = quiet_nan()
         return
      end if
      if (abs(flo) <= tiny(1.0_dp)) then
         root = lo
         return
      else if (abs(fhi) <= tiny(1.0_dp)) then
         root = hi
         return
      end if

      do i = 1, iterations
         mid = 0.5_dp*(lo+hi)
         fmid = fun%evaluate(mid)
         if (.not. ieee_is_finite(fmid)) then
            root = quiet_nan()
            return
         end if
         if (abs(fmid) <= tiny(1.0_dp) .or. abs(hi-lo) <= atol+rtol*abs(mid)) exit
         if (flo*fmid <= 0.0_dp) then
            hi = mid
            fhi = fmid
         else
            lo = mid
            flo = fmid
         end if
      end do
      root = 0.5_dp*(lo+hi)
   end function solve_bisection

end module risk_numerics
