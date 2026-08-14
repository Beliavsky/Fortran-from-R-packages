! SPDX-License-Identifier: GPL-2.0-only
program test_reproducibility
   use mco, only : dp, nsga2_options, nsga2_result, nsga2_optimize
   implicit none
   type(nsga2_options)::opt
   type(nsga2_result)::a,b
   opt%population_size=40; opt%generations=20; opt%seed=2026
   call nsga2_optimize(obj,2,2,[0.0_dp,0.0_dp],[1.0_dp,1.0_dp],a,opt)
   call nsga2_optimize(obj,2,2,[0.0_dp,0.0_dp],[1.0_dp,1.0_dp],b,opt)
   if(maxval(abs(a%par-b%par))>0.0_dp) error stop "seed not reproducible"
   if(any(a%rank/=b%rank)) error stop "rank not reproducible"
   print '(a)', 'test_reproducibility: PASS'
contains
   subroutine obj(x,f)
      real(dp),intent(in)::x(:); real(dp),intent(out)::f(:)
      f=[sum(x*x),sum((x-1.0_dp)**2)]
   end subroutine
end program
