! SPDX-License-Identifier: GPL-2.0-or-later
module test_grad_functions
   use numderiv, only : dp
   implicit none
contains
   function sum_sin(z) result(value)
      real(dp), intent(in) :: z(:)
      real(dp) :: value
      value = sum(sin(z))
   end function sum_sin

   function sum_sin_complex(z) result(value)
      complex(dp), intent(in) :: z(:)
      complex(dp) :: value
      value = sum(sin(z))
   end function sum_sin_complex

   function sin_vector(z) result(value)
      real(dp), intent(in) :: z(:)
      real(dp), allocatable :: value(:)
      value = sin(z)
   end function sin_vector

   function sin_vector_complex(z) result(value)
      complex(dp), intent(in) :: z(:)
      complex(dp), allocatable :: value(:)
      value = sin(z)
   end function sin_vector_complex

   function boundary_function(z) result(value)
      real(dp), intent(in) :: z(:)
      real(dp) :: value
      value = 1.0_dp - (1.0_dp - z(1)) + sin(z(2))
   end function boundary_function

end module test_grad_functions

program test_grad
   use test_grad_functions
   use numderiv, only : dp, deriv_options, nd_success, nd_invalid_argument, grad, grad_complex, &
      grad_elementwise, grad_elementwise_complex
   implicit none

   real(dp) :: x(5), g(5), exact(5), xb(2), gb(2)
   type(deriv_options) :: opts
   integer :: status

   x = [0.0_dp, 0.2_dp, 0.7_dp, 1.3_dp, 2.0_dp]
   exact = cos(x)

   call grad(sum_sin, x, g, status=status)
   call check(status == nd_success, 'Richardson gradient status')
   call check_close(g, exact, 2.0e-10_dp, 'Richardson gradient')

   call grad(sum_sin, x, g, method='simple', status=status)
   call check_close(g, exact, 6.0e-5_dp, 'simple gradient')

   call grad_complex(sum_sin_complex, x, g, status=status)
   call check_close(g, exact, 2.0e-15_dp, 'complex gradient')

   call grad_elementwise(sin_vector, x, g, status=status)
   call check_close(g, exact, 2.0e-10_dp, 'elementwise Richardson gradient')

   call grad_elementwise_complex(sin_vector_complex, x, g, status=status)
   call check_close(g, exact, 2.0e-15_dp, 'elementwise complex gradient')

   xb = [1.0_dp, 0.3_dp]
   call grad(boundary_function, xb, gb, side=[-1, 0], status=status)
   call check(status == nd_success, 'one-sided gradient status')
   call check(abs(gb(1) - 1.0_dp) < 2.0e-8_dp, 'one-sided first derivative')
   call check(abs(gb(2) - cos(xb(2))) < 2.0e-10_dp, 'central second derivative')

   opts = deriv_options(r=1)
   call grad(sum_sin, x, g, options=opts, status=status)
   call check(status == nd_success, 'r=1 status')
   call check_close(g, exact, 2.0e-8_dp, 'r=1 central difference')

   opts = deriv_options(eps=-1.0_dp)
   call grad(sum_sin, x, g, options=opts, status=status)
   call check(status == nd_invalid_argument, 'invalid option status')

   print '(a)', 'test_grad: PASS'

contains

   subroutine check(condition, label)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: label
      if (.not. condition) then
         write(*, '(a,1x,a)') 'FAIL:', trim(label)
         error stop 1
      end if
   end subroutine check

   subroutine check_close(actual, expected, tolerance, label)
      real(dp), intent(in) :: actual(:), expected(:), tolerance
      character(len=*), intent(in) :: label
      call check(maxval(abs(actual - expected)) <= tolerance, label)
   end subroutine check_close

end program test_grad
