! SPDX-License-Identifier: Artistic-2.0
program test_varma_varx
   use mts
   use test_support
   implicit none
   integer,parameter::n=700,k=2
   real(dp)::intercept(k),phi(k,k,1),theta(k,k,1),sigma(k,k),x(n,k),innov(n,k)
   real(dp)::exog(n,1),z(n,k),noise(k),true_beta(k,1),true_phi(k,k)
   real(dp),allocatable::res(:,:),psi(:,:,:),covs(:,:,:),fc(:,:),fcov(:,:,:)
   real(dp),allocatable::coef(:,:),mlres(:,:),mlsigma(:,:),se(:,:),aic(:,:),bic(:,:),hq(:,:)
   type(varma_model)::model
   type(varx_model)::vx
   integer::i,st,bp,bm,istat

   call set_random_seed(12345)
   intercept=[0.03_dp,-0.02_dp]
   phi(:,:,1)=reshape([0.35_dp,0.04_dp,-0.06_dp,0.25_dp],[k,k])
   theta(:,:,1)=reshape([0.20_dp,0.00_dp,0.02_dp,0.15_dp],[k,k])
   sigma=reshape([0.10_dp,0.02_dp,0.02_dp,0.08_dp],[k,k])
   call simulate_varma(n,intercept,phi,theta,sigma,x,innov,burn_in=300)
   allocate(res(n,k));call varma_residuals(x,intercept,phi,theta,res,st)
   call assert_true(st==2,'VARMA residual start')
   call assert_true(sum(res(st:n,:)**2)>0.0_dp,'VARMA residual values')
   call varma_psi_weights(phi,theta,4,psi)
   call assert_matrix_close(psi(:,:,1),phi(:,:,1)-theta(:,:,1),1.0e-12_dp,'VARMA psi first lag')
   call varma_covariance(phi,theta,sigma,3,covs,truncation=200)
   call assert_true(all(shape(covs)==[k,k,4]),'VARMA covariance shape')

   call fit_varma(x(1:450,:),1,1,model,max_iterations=80,tolerance=1.0e-6_dp)
   call assert_true(allocated(model%sigma).and.model%iterations>=0,'VARMA fit output')
   call assert_finite([model%log_likelihood,model%aic,model%bic],'finite VARMA criteria')
   call predict_varma(model,x(1:450,:),3,fc,fcov)
   call assert_true(all(shape(fc)==[3,k]),'VARMA forecast shape')

   true_phi=reshape([0.40_dp,0.05_dp,-0.08_dp,0.20_dp],[k,k])
   true_beta(:,1)=[0.7_dp,-0.4_dp]
   exog(1,1)=0.0_dp;z(1,:)=[0.1_dp,-0.1_dp]
   do i=2,n
      exog(i,1)=0.7_dp*exog(i-1,1)+0.2_dp*sin(0.1_dp*real(i,dp))
      call random_multivariate_normal([0.0_dp,0.0_dp],0.01_dp*eye(k),noise)
      z(i,:)=[0.05_dp,-0.03_dp]+matmul(true_phi,z(i-1,:))+matmul(true_beta,exog(i,:))+noise
   end do
   call fit_varx(z,1,exog,0,vx)
   call assert_true(vx%status==mts_success,'VARX fit status')
   call assert_true(maxval(abs(vx%phi(:,:,1)-true_phi))<0.08_dp,'VARX AR recovery')
   call assert_true(maxval(abs(vx%beta(:,:,0)-true_beta))<0.08_dp,'VARX exogenous recovery')
   call select_varx_order(z(1:300,:),exog(1:300,:),2,1,aic,bic,hq,bp,bm)
   call assert_true(bp>=0.and.bp<=2.and.bm>=0.and.bm<=1,'VARX order selection')

   call multivariate_linear_model(z(2:n,:),exog(2:n,:),coef,mlres,mlsigma,se,status=istat)
   call assert_true(istat==mts_success.and.all(shape(coef)==[2,k]),'multivariate regression')

   print '(a)','test_varma_varx: PASS'
end program test_varma_varx
