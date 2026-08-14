! SPDX-License-Identifier: GPL-2.0-only
program example_binh1
   use mco, only : dp, nsga2_options, nsga2_result, nsga2_optimize, binh1
   implicit none
   type(nsga2_options) :: opt
   type(nsga2_result) :: res
   opt%population_size=80; opt%generations=100; opt%seed=123
   call nsga2_optimize(binh1,2,2,[-5.0_dp,-5.0_dp],[10.0_dp,10.0_dp],res,opt)
   print '(a,i0)', 'Pareto points: ',count(res%pareto_optimal)
   print '(a,2f12.6)', 'Objective ranges: ',minval(res%value(1,:)),maxval(res%value(1,:))
end program
