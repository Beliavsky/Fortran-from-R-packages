! SPDX-License-Identifier: Artistic-2.0
program demo_mts
   use mts
   implicit none
   integer,parameter::n=500,k=2
   real(dp)::x(n,k),intercept(k),phi(k,k,1),sigma(k,k)
   real(dp),allocatable::forecast(:,:),fcov(:,:,:)
   type(var_model)::model

   call set_random_seed(20260801)
   intercept=[0.05_dp,-0.02_dp]
   phi(:,:,1)=reshape([0.55_dp,0.08_dp,-0.12_dp,0.30_dp],[k,k])
   sigma=reshape([0.20_dp,0.04_dp,0.04_dp,0.12_dp],[k,k])
   call simulate_var(n,intercept,phi,sigma,x)
   call fit_var(x,1,model)
   call predict_var(model,x,4,forecast,fcov)

   print '(a)','MTS modern Fortran demonstration'
   print '(a,i0)','status: ',model%status
   print '(a,2(1x,f9.5))','estimated intercept:',model%intercept
   print '(a)','# estimated lag-one matrix'
   print '(2(1x,f9.5))',model%phi(:,:,1)
   print '(a)','# four-step forecasts'
   print '(2(1x,f9.5))',transpose(forecast)
end program demo_mts
