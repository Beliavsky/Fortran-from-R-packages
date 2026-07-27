! SPDX-License-Identifier: GPL-3.0-only
program fhmm_demo
   use fhmm
   implicit none
   type(hmm_parameters) :: par
   type(hmm_simulation) :: sim
   type(inference_result) :: inf
   integer, allocatable :: states(:)
   par%distribution=dist_normal
   allocate(par%gamma(2,2),par%mu(2),par%sigma(2),par%df(2))
   par%gamma=reshape([0.97_dp,0.06_dp,0.03_dp,0.94_dp],[2,2])
   par%mu=[-1.0_dp,1.0_dp];par%sigma=[0.7_dp,0.9_dp];par%df=10.0_dp
   sim=simulate_hmm_model(par,250,seed=2026)
   inf=forward_backward(sim%observations,par)
   states=viterbi_decode(sim%observations,par)
   print '(a,f12.4)','log likelihood: ',inf%log_likelihood
   print '(a,2f10.4)','state frequencies: ',real(count(states==1),dp)/250.0_dp,real(count(states==2),dp)/250.0_dp
end program fhmm_demo
