! SPDX-License-Identifier: GPL-2.0-only
program demo_mco
   use mco, only : dp, nsga2_options, nsga2_result, nsga2_optimize, fonseca2, &
      dominated_hypervolume
   implicit none
   type(nsga2_options) :: opt
   type(nsga2_result) :: res
   opt%population_size=100; opt%generations=150; opt%seed=2026
   call nsga2_optimize(fonseca2,2,2,[-4.0_dp,-4.0_dp],[4.0_dp,4.0_dp],res,opt)
   print '(a,i0)', 'Evaluations: ',res%evaluations
   print '(a,i0)', 'Pareto points: ',count(res%pareto_optimal)
   print '(a,f10.6)', 'Dominated hypervolume to (1,1): ',dominated_hypervolume(res%value,[1.0_dp,1.0_dp])
end program
