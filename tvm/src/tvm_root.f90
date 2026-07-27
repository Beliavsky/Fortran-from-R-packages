! SPDX-License-Identifier: MIT
! Copyright (c) 2014 Juan Manuel Truppia
! Modern Fortran translation of the tvm package.
module tvm_root
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use tvm_kinds, only : dp
   implicit none
   private

   integer, parameter, public :: root_success = 0
   integer, parameter, public :: root_not_bracketed = 1
   integer, parameter, public :: root_max_iterations = 2
   integer, parameter, public :: root_nonfinite = 3

   abstract interface
      function scalar_function(x) result(y)
         import dp
         real(dp), intent(in) :: x
         real(dp) :: y
      end function scalar_function
   end interface

   public :: bisect_root

contains

   subroutine bisect_root(f, lower, upper, root, status, tol, max_iterations)
      procedure(scalar_function) :: f
      real(dp), intent(in) :: lower, upper
      real(dp), intent(out) :: root
      integer, intent(out) :: status
      real(dp), intent(in), optional :: tol
      integer, intent(in), optional :: max_iterations
      real(dp) :: a, b, c, fa, fb, fc, xtol
      integer :: iter, nmax

      xtol = 1.0e-10_dp
      if (present(tol)) xtol = max(tol, epsilon(1.0_dp))
      nmax = 200
      if (present(max_iterations)) nmax = max(1, max_iterations)

      a = lower
      b = upper
      if (b <= a) then
         root = a
         status = root_not_bracketed
         return
      end if

      fa = f(a)
      fb = f(b)
      if (.not. finite_value(fa) .or. .not. finite_value(fb)) then
         root = 0.5_dp * (a + b)
         status = root_nonfinite
         return
      end if
      if (is_zero(fa)) then
         root = a
         status = root_success
         return
      end if
      if (is_zero(fb)) then
         root = b
         status = root_success
         return
      end if
      if (same_sign(fa, fb)) then
         root = 0.5_dp * (a + b)
         status = root_not_bracketed
         return
      end if

      do iter = 1, nmax
         c = 0.5_dp * (a + b)
         fc = f(c)
         if (.not. finite_value(fc)) then
            root = c
            status = root_nonfinite
            return
         end if
         if (is_zero(fc) .or. abs(b - a) <= 2.0_dp * xtol * max(1.0_dp, abs(c))) then
            root = c
            status = root_success
            return
         end if
         if (same_sign(fa, fc)) then
            a = c
            fa = fc
         else
            b = c
            fb = fc
         end if
      end do

      root = 0.5_dp * (a + b)
      status = root_max_iterations
   end subroutine bisect_root

   pure logical function same_sign(a, b) result(answer)
      real(dp), intent(in) :: a, b
      answer = (a > 0.0_dp .and. b > 0.0_dp) .or. (a < 0.0_dp .and. b < 0.0_dp)
   end function same_sign

   pure logical function is_zero(x) result(answer)
      real(dp), intent(in) :: x
      answer = abs(x) <= tiny(1.0_dp)
   end function is_zero

   pure logical function finite_value(x) result(answer)
      real(dp), intent(in) :: x
      answer = ieee_is_finite(x)
   end function finite_value

end module tvm_root
