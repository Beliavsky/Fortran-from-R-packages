! SPDX-License-Identifier: Artistic-2.0
program test_diagnostics_factor
   use mts
   use test_support
   implicit none
   integer,parameter::n=600,k=4
   real(dp)::x(n,k),factor(n),eps(k),h(k,2),y(n)
   real(dp),allocatable::ccm(:,:,:),pv(:),omega(:,:),factors(:,:),psi(:,:),loadings(:,:)
   real(dp),allocatable::coef(:),forecast(:),swload(:,:),swfactor(:,:),priorp(:,:),priors(:,:)
   real(dp)::mse
   type(diagnostic_result)::port,arch,rankarch
   type(factor_result)::pca,apca
   type(bvar_result)::bvar
   integer::i,j,istat,d

   call set_random_seed(4567)
   do i=1,n
      call random_multivariate_normal([0.0_dp,0.0_dp,0.0_dp,0.0_dp],eye(k),eps)
      factor(i)=sin(0.035_dp*real(i,dp))+0.25_dp*eps(1)
      x(i,:)=[1.0_dp,0.8_dp,-0.6_dp,0.4_dp]*factor(i)+0.20_dp*eps
      y(i)=0.5_dp+1.2_dp*factor(i)+0.1_dp*eps(2)
   end do
   call cross_correlation_matrices(x,5,ccm,pv,istat)
   call assert_true(istat==mts_success.and.all(shape(ccm)==[k,k,6]),'cross-correlation matrices')
   call multivariate_portmanteau(x,5,port,adjusted=.true.)
   call multivariate_arch_test(x,3,arch)
   call rank_arch_test(x,3,rankarch)
   call assert_true(port%status==mts_success.and.arch%status==mts_success.and.rankarch%status==mts_success,'diagnostic tests')
   call assert_true(port%p_value>=0.0_dp.and.port%p_value<=1.0_dp,'portmanteau p-value')

   call principal_components(x,1,pca,standardize=.true.)
   call assert_true(pca%status==mts_success.and.all(shape(pca%loadings)==[k,1]),'PCA')
   call assert_true(pca%eigenvalues(1)>pca%eigenvalues(2),'PCA eigen ordering')
   call asymptotic_pca(transpose(x(1:4,1:4)),1,apca)
   call assert_true(apca%status==mts_success,'asymptotic PCA')

   h=reshape([1.0_dp,0.0_dp,1.0_dp,0.0_dp,0.0_dp,1.0_dp,0.0_dp,1.0_dp],[k,2])
   call constrained_factor_model(x,h,1,omega,factors,psi,loadings,istat)
   call assert_true(istat==mts_success.and.all(shape(loadings)==[k,1]),'constrained factor model')
   call stock_watson_forecast(y,x,450,1,coef,forecast,mse,swload,swfactor,istat)
   call assert_true(istat==mts_success.and.size(forecast)==n-450.and.mse>=0.0_dp,'Stock-Watson forecast')

   d=1+k
   priorp=0.1_dp*eye(d);priors=eye(k)
   call fit_bvar(x,1,priorp,priors,k+2,bvar)
   call assert_true(bvar%status==mts_success.and.all(shape(bvar%posterior_mean)==[d,k]),'BVAR')

   print '(a)','test_diagnostics_factor: PASS'
end program test_diagnostics_factor
