! SPDX-License-Identifier: GPL-3.0-or-later
program test_pin_estimation
   use pinstimation
   implicit none
   type(pin_parameters) :: truth
   type(pin_result) :: fit
   type(bayes_pin_result) :: bayes
   type(trade_counts) :: data
   integer, allocatable :: states(:)
   real(dp) :: true_pin

   truth = pin_parameters(0.38_dp,0.42_dp,16.0_dp,20.0_dp,18.0_dp)
   call simulate_pin(450,truth,data,states,seed=8127)
   call fit_pin(data,fit,max_iterations=1800,tolerance=1.0e-7_dp)
   if (.not. finite_number(fit%log_likelihood)) error stop 'nonfinite PIN fit'
   true_pin = pin_value(truth)
   if (abs(fit%pin-true_pin)>0.09_dp) then
      write(*,*) 'fit PIN too far from truth',fit%pin,true_pin
      error stop 1
   end if
   if (minval(fit%posteriors)<-1.0e-12_dp .or. maxval(fit%posteriors)>1.0_dp+1.0e-12_dp) error stop 'invalid posterior'
   call fit_pin_bayes(data,bayes,sweeps=500,burnin=150,thin=2,initial=fit%parameters,proposal_scale=0.08_dp,seed=77)
   if (size(bayes%draws,1)/=175) error stop 'wrong number of Bayesian draws'
   if (bayes%acceptance_rate<=0.0_dp .or. bayes%acceptance_rate>=1.0_dp) error stop 'invalid acceptance rate'
   print '(a)', 'test_pin_estimation: PASS'
end program test_pin_estimation
