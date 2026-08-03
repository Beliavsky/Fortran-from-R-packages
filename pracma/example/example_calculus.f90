! SPDX-License-Identifier: GPL-3.0-or-later
module example_calculus_functions
   use pracma_kinds, only : dp
   implicit none
contains
   function integrand(x) result(y)
      real(dp),intent(in)::x
      real(dp)::y
      y=exp(-x*x)
   end function integrand

   function objective(x) result(y)
      real(dp),intent(in)::x(:)
      real(dp)::y
      y=(1.0_dp-x(1))**2+100.0_dp*(x(2)-x(1)*x(1))**2
   end function objective
end module example_calculus_functions

program example_calculus
   use pracma
   use example_calculus_functions
   implicit none
   type(quadrature_result) :: q
   type(optimization_result) :: opt

   q=integral(integrand,0.0_dp,1.0_dp)
   opt=fminsearch(objective,[-1.2_dp,1.0_dp],max_iter=5000)
   print '(a,f14.10)','integral exp(-x^2), 0..1: ',q%value
   print '(a,2f12.6)','Rosenbrock minimizer: ',opt%x
   print '(a,es12.4)','objective: ',opt%value
end program example_calculus
