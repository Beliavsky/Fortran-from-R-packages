! SPDX-License-Identifier: GPL-3.0-or-later
program adjpin_estimation
   use pinstimation
   implicit none
   type(adjpin_parameters) :: truth
   type(adjpin_restrictions) :: restrictions
   type(adjpin_result) :: estimate
   type(trade_counts) :: data
   integer, allocatable :: state(:)

   truth=adjpin_parameters(0.30_dp,0.40_dp,0.20_dp,0.20_dp,18.0_dp,18.0_dp,12.0_dp,12.0_dp,5.0_dp,5.0_dp)
   restrictions%equal_theta=.true.
   restrictions%equal_eps=.true.
   restrictions%equal_mu=.true.
   restrictions%equal_d=.true.
   call simulate_adjpin(300,truth,data,state,seed=910)
   call fit_adjpin(data,estimate,restrictions=restrictions,method='ML',max_iterations=1800)
   print '(a,f12.6)', 'AdjPIN: ',estimate%adjpin
   print '(a,f12.6)', 'PSOS:   ',estimate%psos
   print '(a,f14.4)', 'log likelihood: ',estimate%log_likelihood
end program adjpin_estimation
