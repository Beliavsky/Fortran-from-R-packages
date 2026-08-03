! SPDX-License-Identifier: Artistic-2.0
program varma_analysis
   use mts
   implicit none
   integer,parameter::n=500,k=2
   real(dp)::x(n,k),mu(k),phi(k,k,1),theta(k,k,1),sigma(k,k)
   real(dp),allocatable::forecast(:,:),covariance(:,:,:)
   type(varma_model)::fit

   call set_random_seed(2202)
   mu=[0.01_dp,-0.01_dp]
   phi(:,:,1)=reshape([0.35_dp,0.02_dp,-0.05_dp,0.25_dp],[k,k])
   theta(:,:,1)=reshape([0.18_dp,0.00_dp,0.03_dp,0.12_dp],[k,k])
   sigma=reshape([0.10_dp,0.02_dp,0.02_dp,0.08_dp],[k,k])
   call simulate_varma(n,mu,phi,theta,sigma,x)
   call fit_varma(x,1,1,fit,max_iterations=100)
   call predict_varma(fit,x,3,forecast,covariance)

   print '(a,i0)','VARMA fit status: ',fit%status
   print '(a,f12.4)','log likelihood: ',fit%log_likelihood
   print '(a,2(1x,f9.5))','one-step forecast:',forecast(1,:)
end program varma_analysis
