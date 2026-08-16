! SPDX-License-Identifier: GPL-2.0-or-later
module gradient_example_functions
   use numderiv, only : dp
   implicit none
contains
   function objective(x) result(value)
      real(dp), intent(in) :: x(:)
      real(dp) :: value
      value = sum((exp(x) - x) * [1.0_dp, 2.0_dp, 3.0_dp]) / 3.0_dp
   end function objective

   function objective_complex(x) result(value)
      complex(dp), intent(in) :: x(:)
      complex(dp) :: value
      value = sum((exp(x) - x) * [1.0_dp, 2.0_dp, 3.0_dp]) / 3.0_dp
   end function objective_complex
end module gradient_example_functions

program gradient_example
   use gradient_example_functions
   use numderiv, only : dp, grad, grad_complex
   implicit none
   real(dp) :: x(3), g_richardson(3), g_complex(3)

   x = [0.1_dp, -0.2_dp, 0.5_dp]
   call grad(objective, x, g_richardson)
   call grad_complex(objective_complex, x, g_complex)

   print '(a,3(1x,es16.8))', 'Richardson:', g_richardson
   print '(a,3(1x,es16.8))', 'complex step:', g_complex
end program gradient_example
