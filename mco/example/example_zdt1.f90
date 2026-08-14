! SPDX-License-Identifier: GPL-2.0-only
program example_zdt1
   use mco, only : dp, nsga2_options, nsga2_result, nsga2_optimize, zdt1
   implicit none
   type(nsga2_options) :: opt
   type(nsga2_result) :: res
   real(dp) :: lower(30),upper(30),err
   lower=0.0_dp; upper=1.0_dp
   opt%population_size=100; opt%generations=300; opt%seed=42
   call nsga2_optimize(zdt1,30,2,lower,upper,res,opt)
   err=maxval(res%value(2,:)-(1.0_dp-sqrt(res%value(1,:))))
   print '(a,es12.4)', 'Maximum absolute vertical gap from ZDT1 front: ',err
end program
