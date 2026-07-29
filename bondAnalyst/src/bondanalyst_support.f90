! bondAnalyst modern Fortran port
! Copyright (C) 2022 MaheshP Kumar
! Copyright (C) 2026 Fortran port contributors
! SPDX-License-Identifier: GPL-3.0-only
module bondanalyst_support
   use, intrinsic :: ieee_arithmetic, only : ieee_quiet_nan, ieee_value, ieee_is_finite
   use bondanalyst_kinds, only : dp
   implicit none
   private

   integer, parameter, public :: ba_success = 0
   integer, parameter, public :: ba_invalid_argument = 1
   integer, parameter, public :: ba_size_mismatch = 2
   integer, parameter, public :: ba_no_root = 3
   integer, parameter, public :: ba_out_of_range = 4

   public :: round_decimal, quiet_nan, set_status, solve_positive_root
   public :: all_finite, valid_discount_base

   abstract interface
      function scalar_function(x) result(value)
         import :: dp
         real(dp), intent(in) :: x
         real(dp) :: value
      end function scalar_function
   end interface

contains

   pure function quiet_nan() result(value)
      real(dp) :: value

      value = ieee_value(0.0_dp, ieee_quiet_nan)
   end function quiet_nan


   pure function round_decimal(x, digits) result(value)
      real(dp), intent(in) :: x
      integer, intent(in) :: digits
      real(dp) :: value
      real(dp) :: scale

      if (.not. ieee_is_finite(x)) then
         value = x
         return
      end if
      scale = 10.0_dp**digits
      value = anint(x*scale)/scale
   end function round_decimal


   pure subroutine set_status(istat, code)
      integer, optional, intent(out) :: istat
      integer, intent(in) :: code

      if (present(istat)) istat = code
   end subroutine set_status


   pure function all_finite(x) result(ok)
      real(dp), intent(in) :: x(:)
      logical :: ok

      ok = all(ieee_is_finite(x))
   end function all_finite


   pure function valid_discount_base(rate) result(ok)
      real(dp), intent(in) :: rate
      logical :: ok

      ok = ieee_is_finite(rate) .and. rate > -1.0_dp
   end function valid_discount_base


   function solve_positive_root(fun, root, tolerance, max_iterations) result(status)
      procedure(scalar_function) :: fun
      real(dp), intent(out) :: root
      real(dp), optional, intent(in) :: tolerance
      integer, optional, intent(in) :: max_iterations
      integer :: status
      integer :: iter, maxit
      real(dp) :: a, b, c, fa, fb, fc, tol

      tol = 1.0e-12_dp
      if (present(tolerance)) tol = max(tolerance, epsilon(1.0_dp))
      maxit = 256
      if (present(max_iterations)) maxit = max(1, max_iterations)

      a = 0.0_dp
      fa = fun(a)
      if (.not. ieee_is_finite(fa)) then
         root = quiet_nan()
         status = ba_invalid_argument
         return
      end if
      if (abs(fa) <= tol) then
         root = 0.0_dp
         status = ba_success
         return
      end if
      if (fa < 0.0_dp) then
         root = quiet_nan()
         status = ba_no_root
         return
      end if

      b = 0.01_dp
      fb = fun(b)
      do while (ieee_is_finite(fb) .and. fb > 0.0_dp .and. b < 1.0e8_dp)
         b = 2.0_dp*b + 0.01_dp
         fb = fun(b)
      end do
      if (.not. ieee_is_finite(fb) .or. fb > 0.0_dp) then
         root = quiet_nan()
         status = ba_no_root
         return
      end if

      do iter = 1, maxit
         c = 0.5_dp*(a+b)
         fc = fun(c)
         if (.not. ieee_is_finite(fc)) then
            root = quiet_nan()
            status = ba_no_root
            return
         end if
         if (abs(fc) <= tol .or. 0.5_dp*(b-a) <= tol*max(1.0_dp, abs(c))) then
            root = c
            status = ba_success
            return
         end if
         if (fc > 0.0_dp) then
            a = c
            fa = fc
         else
            b = c
            fb = fc
         end if
      end do

      root = 0.5_dp*(a+b)
      status = ba_no_root
   end function solve_positive_root

end module bondanalyst_support
