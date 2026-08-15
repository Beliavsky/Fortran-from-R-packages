! rmutil computational translation
! Copyright (C) 1998-2001 J.K. Lindsey
! Copyright (C) 2026 OpenAI (modern Fortran translation)
! SPDX-License-Identifier: GPL-2.0-or-later
module rmutil_ode
   use rmutil_kinds, only : dp
   implicit none
   private
   public :: runge_kutta
   abstract interface
      function ode_scalar(y,x) result(v)
         import dp
         real(dp), intent(in) :: y, x
         real(dp) :: v
      end function ode_scalar
   end interface
contains
   function runge_kutta(f, initial, x) result(y)
      procedure(ode_scalar) :: f
      real(dp), intent(in) :: initial, x(:)
      real(dp), allocatable :: y(:)
      real(dp) :: h, f1, f2, f3, f4
      integer :: i
      allocate(y(size(x)))
      y(1) = initial
      do i = 1, size(x)-1
         h = x(i+1)-x(i)
         f1 = h*f(y(i),x(i))
         f2 = h*f(y(i)+f1/2.0_dp,x(i)+h/2.0_dp)
         f3 = h*f(y(i)+f2/2.0_dp,x(i)+h/2.0_dp)
         f4 = h*f(y(i)+f3,x(i)+h)
         y(i+1) = y(i) + (f1+2.0_dp*f2+2.0_dp*f3+f4)/6.0_dp
      end do
   end function runge_kutta
end module rmutil_ode
