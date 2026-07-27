! SPDX-License-Identifier: GPL-3.0-only
program test_estimation
   use fhmm
   implicit none
   type(hmm_parameters) :: truth,start
   type(hmm_simulation) :: sim
   type(hmm_fit_result) :: fit
   type(fit_options) :: opt

   truth%distribution=dist_normal
   allocate(truth%gamma(2,2),truth%mu(2),truth%sigma(2),truth%df(2))
   truth%gamma=reshape([0.96_dp,0.05_dp,0.04_dp,0.95_dp],[2,2])
   truth%mu=[-1.5_dp,1.5_dp];truth%sigma=[0.6_dp,0.7_dp];truth%df=10.0_dp
   sim=simulate_hmm_model(truth,350,seed=918)
   start=truth
   start%mu=start%mu+[0.2_dp,-0.2_dp]
   opt%runs=1;opt%max_iterations=600;opt%x_tolerance=2.0e-5_dp;opt%f_tolerance=2.0e-5_dp
   opt%compute_hessian=.false.;opt%seed=11
   fit=fit_hmm(sim%observations,2,dist_normal,opt,start)
   if(.not.fit%ok)then;write(*,*)trim(fit%message);error stop 1;end if
   if(fit%log_likelihood<hmm_log_likelihood(sim%observations,start)-1.0e-6_dp)error stop 1
   if(maxval(abs(sum(fit%parameters%gamma,dim=2)-1.0_dp))>1.0e-12_dp)error stop 1
   if(any(fit%parameters%sigma<=0.0_dp))error stop 1
   if(size(fit%decoding)/=350)error stop 1
   print '(a)','test_estimation: PASS'
end program test_estimation
