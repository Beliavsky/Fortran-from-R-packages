! SPDX-License-Identifier: GPL-2.0-or-later
module test_jacobian_functions
   use numderiv, only : dp
   implicit none
contains
   function vector_function(z) result(value)
      real(dp), intent(in) :: z(:)
      real(dp), allocatable :: value(:)
      allocate(value(3))
      value = [z(1) ** 2 + z(2), sin(z(1) * z(2)), exp(z(2))]
   end function vector_function

   function vector_function_complex(z) result(value)
      complex(dp), intent(in) :: z(:)
      complex(dp), allocatable :: value(:)
      allocate(value(3))
      value = [z(1) ** 2 + z(2), sin(z(1) * z(2)), exp(z(2))]
   end function vector_function_complex

end module test_jacobian_functions

program test_jacobian
   use test_jacobian_functions
   use numderiv, only : dp, nd_success, jacobian, jacobian_complex
   implicit none

   real(dp) :: x(2), exact(3, 2)
   real(dp), allocatable :: jac(:, :)
   integer :: status

   x = [0.7_dp, -0.4_dp]
   exact(1, :) = [2.0_dp * x(1), 1.0_dp]
   exact(2, :) = [x(2) * cos(x(1) * x(2)), x(1) * cos(x(1) * x(2))]
   exact(3, :) = [0.0_dp, exp(x(2))]

   call jacobian(vector_function, x, jac, status=status)
   call check(status == nd_success, 'Richardson Jacobian status')
   call check_close(jac, exact, 2.0e-10_dp, 'Richardson Jacobian')

   call jacobian(vector_function, x, jac, method='simple', side=[1, -1], status=status)
   call check(status == nd_success, 'simple Jacobian status')
   call check_close(jac, exact, 1.2e-4_dp, 'simple Jacobian')

   call jacobian_complex(vector_function_complex, x, jac, status=status)
   call check(status == nd_success, 'complex Jacobian status')
   call check_close(jac, exact, 2.0e-15_dp, 'complex Jacobian')

   print '(a)', 'test_jacobian: PASS'

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
      real(dp), intent(in) :: actual(:, :), expected(:, :), tolerance
      character(len=*), intent(in) :: label
      call check(maxval(abs(actual - expected)) <= tolerance, label)
   end subroutine check_close

end program test_jacobian
