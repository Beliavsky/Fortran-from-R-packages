! SPDX-License-Identifier: GPL-3.0-or-later
program mpin_estimation
   use pinstimation
   implicit none
   type(mpin_parameters) :: truth
   type(mpin_result) :: estimate
   type(trade_counts) :: data
   integer, allocatable :: state(:)

   allocate(truth%alpha(2),truth%delta(2),truth%mu(2))
   truth%alpha=[0.18_dp,0.22_dp]
   truth%delta=[0.30_dp,0.65_dp]
   truth%mu=[10.0_dp,25.0_dp]
   truth%eps_b=18.0_dp
   truth%eps_s=19.0_dp
   call simulate_mpin(350,truth,data,state,seed=55)
   call fit_mpin_ecm(data,2,estimate,max_iterations=60)
   print '(a,2f12.5)', 'alpha: ',estimate%parameters%alpha
   print '(a,2f12.5)', 'delta: ',estimate%parameters%delta
   print '(a,2f12.5)', 'mu:    ',estimate%parameters%mu
   print '(a,f12.6)', 'MPIN:  ',estimate%mpin
end program mpin_estimation
