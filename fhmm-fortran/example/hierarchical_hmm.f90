! SPDX-License-Identifier: GPL-3.0-only
program hierarchical_hmm
   use fhmm
   implicit none
   type(hhmm_parameters) :: par
   type(hhmm_simulation) :: sim
   type(hhmm_inference_result) :: inf
   integer :: s
   par%coarse%distribution=dist_normal
   allocate(par%coarse%gamma(2,2),par%coarse%mu(2),par%coarse%sigma(2),par%coarse%df(2),par%fine(2))
   par%coarse%gamma=reshape([0.95_dp,0.08_dp,0.05_dp,0.92_dp],[2,2])
   par%coarse%mu=[-1.0_dp,1.0_dp];par%coarse%sigma=0.6_dp;par%coarse%df=10.0_dp
   do s=1,2
      par%fine(s)%distribution=dist_normal
      allocate(par%fine(s)%gamma(2,2),par%fine(s)%mu(2),par%fine(s)%sigma(2),par%fine(s)%df(2))
      par%fine(s)%gamma=reshape([0.9_dp,0.1_dp,0.1_dp,0.9_dp],[2,2])
      par%fine(s)%sigma=0.5_dp;par%fine(s)%df=10.0_dp
   end do
   par%fine(1)%mu=[-1.5_dp,-0.3_dp];par%fine(2)%mu=[0.3_dp,1.5_dp]
   sim=simulate_hhmm_model(par,[5,4,6,5,3],seed=22)
   inf=hhmm_forward_backward(sim%coarse_observations,sim%fine_observations,sim%chunk_lengths,par)
   print '(a,f12.6)','hierarchical log likelihood = ',inf%log_likelihood
end program hierarchical_hmm
