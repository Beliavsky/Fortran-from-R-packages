! SPDX-License-Identifier: GPL-2.0-only
program example_constrained
   use mco, only : dp, nsga2_options, nsga2_result, nsga2_optimize, binh2, binh2_constraints
   implicit none
   type(nsga2_options) :: opt
   type(nsga2_result) :: res
   opt%population_size=80; opt%generations=120; opt%seed=321
   call nsga2_optimize(binh2,2,2,[0.0_dp,0.0_dp],[5.0_dp,3.0_dp],res,opt,binh2_constraints,2)
   print '(a,i0)', 'Feasible Pareto points: ',count(res%pareto_optimal)
   print '(a,es12.4)', 'Maximum survivor violation: ',maxval(res%violation)
end program
