! SPDX-License-Identifier: GPL-3.0-only
program basic_hmm
   use fhmm
   implicit none
   type(hmm_parameters) :: par
   type(inference_result) :: inf
   real(dp) :: y(8)
   par%distribution=dist_student_t
   allocate(par%gamma(2,2),par%mu(2),par%sigma(2),par%df(2))
   par%gamma=reshape([0.95_dp,0.10_dp,0.05_dp,0.90_dp],[2,2])
   par%mu=[-0.5_dp,0.8_dp];par%sigma=[0.7_dp,1.2_dp];par%df=[6.0_dp,8.0_dp]
   y=[-0.8_dp,-0.4_dp,-0.6_dp,0.2_dp,1.1_dp,0.7_dp,1.4_dp,0.4_dp]
   inf=forward_backward(y,par)
   print '(a,f12.6)','log likelihood = ',inf%log_likelihood
   print '(a,2f10.6)','final filtered probabilities = ',inf%filtered(:,size(y))
end program basic_hmm
