! SPDX-License-Identifier: Artistic-2.0
program var_analysis
   use mts
   implicit none
   integer,parameter::n=800,k=3
   real(dp)::x(n,k),mu(k),phi(k,k,1),sigma(k,k)
   real(dp),allocatable::irf(:,:,:),fevd(:,:,:)
   type(var_model)::fit
   type(order_selection_result)::orders

   call set_random_seed(1101)
   mu=[0.02_dp,0.00_dp,-0.01_dp]
   phi(:,:,1)=reshape([0.50_dp,0.03_dp,0.00_dp,-0.10_dp,0.35_dp,0.04_dp,0.05_dp,-0.02_dp,0.25_dp],[k,k])
   sigma=0.08_dp*eye(k);sigma(1,2)=0.01_dp;sigma(2,1)=0.01_dp
   call simulate_var(n,mu,phi,sigma,x)
   call select_var_order(x,4,orders)
   call fit_var(x,orders%bic_order,fit)
   call var_impulse_response(fit%phi,fit%sigma,10,irf)
   call forecast_error_variance_decomposition(fit%phi,fit%sigma,10,fevd)

   print '(a,i0)','BIC-selected VAR order: ',orders%bic_order
   print '(a,3(1x,f9.5))','intercept:',fit%intercept
   print '(a,3(1x,f9.5))','10-step FEVD for variable 1:',fevd(1,:,10)
end program var_analysis
