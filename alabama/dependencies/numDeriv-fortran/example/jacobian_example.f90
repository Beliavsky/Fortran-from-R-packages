! SPDX-License-Identifier: GPL-2.0-or-later
module jacobian_example_functions
   use numderiv, only : dp
   implicit none
contains
   function equations(x) result(value)
      real(dp), intent(in) :: x(:)
      real(dp), allocatable :: value(:)
      value = [x(1) ** 2 + x(2), sin(x(1) * x(2)), exp(x(2))]
   end function equations
end module jacobian_example_functions

program jacobian_example
   use jacobian_example_functions, only : equations
   use numderiv, only : dp, jacobian
   implicit none
   real(dp), allocatable :: jac(:, :)
   real(dp) :: x(2)
   integer :: i

   x = [0.7_dp, -0.4_dp]
   call jacobian(equations, x, jac)
   print '(a)', 'Jacobian:'
   do i = 1, size(jac, 1)
      print '(*(1x,es16.8))', jac(i, :)
   end do
end program jacobian_example
