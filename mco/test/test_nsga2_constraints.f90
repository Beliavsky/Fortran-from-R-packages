! SPDX-License-Identifier: GPL-2.0-only
program test_nsga2_constraints
   use mco, only : dp, nsga2_options, nsga2_result, nsga2_optimize
   implicit none
   type(nsga2_options) :: opt
   type(nsga2_result) :: res
   real(dp) :: lower(2),upper(2), sums(60)
   lower=0.0_dp; upper=1.0_dp
   opt%population_size=60; opt%generations=100; opt%seed=99
   call nsga2_optimize(obj,2,2,lower,upper,res,opt,con,1)
   if(res%status/=0) error stop trim(res%message)
   if(maxval(res%violation)>1.0e-10_dp) error stop "infeasible survivors"
   sums=sum(res%par,dim=1)
   if(minval(sums)<1.0_dp-1.0e-10_dp) error stop "constraint convention"
   if(minval(sums)>1.03_dp) error stop "front did not approach boundary"
   print '(a)', 'test_nsga2_constraints: PASS'
contains
   subroutine obj(x,f)
      real(dp),intent(in)::x(:); real(dp),intent(out)::f(:); f=x
   end subroutine
   subroutine con(x,g)
      real(dp),intent(in)::x(:); real(dp),intent(out)::g(:); g(1)=sum(x)-1.0_dp
   end subroutine
end program
