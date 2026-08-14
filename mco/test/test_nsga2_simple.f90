! SPDX-License-Identifier: GPL-2.0-only
program test_nsga2_simple
   use mco, only : dp, nsga2_options, nsga2_result, nsga2_optimize
   implicit none
   type(nsga2_options) :: opt
   type(nsga2_result) :: res
   real(dp) :: lower(1),upper(1)
   lower=0.0_dp; upper=1.0_dp
   opt%population_size=80; opt%generations=80; opt%seed=4321
   call nsga2_optimize(obj,1,2,lower,upper,res,opt)
   if(res%status/=0) error stop trim(res%message)
   if(count(res%pareto_optimal)<70) error stop "too few Pareto points"
   if(minval(res%par(1,:))>0.08_dp) error stop "left endpoint not explored"
   if(maxval(res%par(1,:))<0.92_dp) error stop "right endpoint not explored"
   if(maxval(abs(res%value(1,:)-res%par(1,:)**2))>1.0e-14_dp) error stop "objective mismatch"
   print '(a)', 'test_nsga2_simple: PASS'
contains
   subroutine obj(x,f)
      real(dp),intent(in)::x(:); real(dp),intent(out)::f(:)
      f=[x(1)**2,(1.0_dp-x(1))**2]
   end subroutine
end program
