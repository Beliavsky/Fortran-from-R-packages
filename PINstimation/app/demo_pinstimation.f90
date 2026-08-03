! SPDX-License-Identifier: GPL-3.0-or-later
program demo_pinstimation
   use pinstimation
   implicit none
   type(pin_parameters) :: truth
   type(pin_result) :: estimate
   type(trade_counts) :: data
   integer,allocatable :: state(:)

   truth=pin_parameters(0.40_dp,0.35_dp,14.0_dp,21.0_dp,19.0_dp)
   call simulate_pin(250,truth,data,state,seed=1)
   call fit_pin(data,estimate,max_iterations=1400)
   print '(a)', 'PINstimation modern Fortran demo'
   print '(a,f10.6)', 'estimated PIN: ',estimate%pin
   print '(a,f10.6)', 'good-news component: ',estimate%pin_good
   print '(a,f10.6)', 'bad-news component:  ',estimate%pin_bad
end program demo_pinstimation
