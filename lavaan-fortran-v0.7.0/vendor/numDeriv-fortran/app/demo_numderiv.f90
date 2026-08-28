! SPDX-License-Identifier: GPL-2.0-or-later
module demo_functions
   use numderiv, only : dp
   implicit none
contains
   function objective(x) result(value)
      real(dp), intent(in) :: x(:)
      real(dp) :: value
      value = (1.0_dp - x(1)) ** 2 + 100.0_dp * (x(2) - x(1) ** 2) ** 2
   end function objective
end module demo_functions

program demo_numderiv
   use demo_functions, only : objective
   use numderiv, only : dp, grad, hessian
   implicit none
   real(dp) :: x(2), gradient(2)
   real(dp), allocatable :: hess(:, :)

   x = [-1.2_dp, 1.0_dp]
   call grad(objective, x, gradient)
   call hessian(objective, x, hess)

   print '(a,2(1x,es15.7))', 'x:', x
   print '(a,2(1x,es15.7))', 'gradient:', gradient
   print '(a)', 'Hessian:'
   print '(2(1x,es15.7))', hess(1, :)
   print '(2(1x,es15.7))', hess(2, :)
end program demo_numderiv
