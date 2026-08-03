! SPDX-License-Identifier: GPL-3.0-or-later
program pin_estimation
   use pinstimation
   implicit none
   type(pin_parameters) :: truth
   type(pin_result) :: estimate
   type(trade_counts) :: data
   integer, allocatable :: state(:)

   truth = pin_parameters(0.35_dp,0.45_dp,15.0_dp,20.0_dp,18.0_dp)
   call simulate_pin(300,truth,data,state,seed=2026)
   call fit_pin(data,estimate,max_iterations=1500)
   print '(a,5f12.5)', 'estimated parameters:', estimate%parameters%alpha, estimate%parameters%delta, &
      estimate%parameters%mu, estimate%parameters%eps_b, estimate%parameters%eps_s
   print '(a,f12.6)', 'PIN: ',estimate%pin
   print '(a,f14.4)', 'log likelihood: ',estimate%log_likelihood
end program pin_estimation
