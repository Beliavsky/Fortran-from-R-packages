! SPDX-License-Identifier: Artistic-2.0
program test_volatility
   use mts
   use test_support
   implicit none
   integer,parameter::n=350,k=2
   real(dp)::r(n,k),z(k),lambda,ll,mean(k),c(k,k),a(k,k),b(k,k)
   real(dp),allocatable::ewma(:,:,:),rho(:,:,:),bcov(:,:,:),bres(:,:),mcov(:,:,:),coef(:,:),comp(:,:)
   real(dp),allocatable::loadings(:,:),eigenvalues(:),transformed(:,:),stdres(:,:),ucor(:,:),gcor(:,:)
   type(dcc_model)::dcc
   type(diagnostic_result)::level,squared
   integer::i,j,istat,it

   call set_random_seed(9090)
   r(1,:)=[0.0_dp,0.0_dp]
   do i=2,n
      call random_multivariate_normal([0.0_dp,0.0_dp],reshape([1.0_dp,0.35_dp,0.35_dp,1.0_dp],[k,k]),z)
      r(i,1)=0.15_dp*r(i-1,1)+z(1)
      r(i,2)=0.10_dp*r(i-1,2)+0.15_dp*r(i-1,1)+z(2)
   end do
   call ewma_covariance(r,0.94_dp,ewma,status=istat)
   call assert_true(istat==mts_success.and.all(shape(ewma)==[k,k,n]),'EWMA covariance')
   do i=1,n
      call assert_true(all([(ewma(j,j,i)>0.0_dp,j=1,k)]),'positive EWMA diagonal')
   end do
   call fit_ewma_lambda(r(1:180,:),lambda,ewma,ll,istat,it)
   call assert_true(lambda>0.0_dp.and.lambda<1.0_dp.and.size(ewma,3)==180,'EWMA lambda fit')
   call assert_finite([ll],'finite EWMA likelihood')

   call dcc_correlations(r,0.04_dp,0.92_dp,rho,status=istat)
   call assert_true(istat==mts_success.and.all(shape(rho)==[k,k,n]),'DCC correlations')
   do i=1,n
      call assert_close(rho(1,1,i),1.0_dp,1.0e-10_dp,'DCC first diagonal')
      call assert_close(rho(2,2,i),1.0_dp,1.0e-10_dp,'DCC second diagonal')
   end do
   call fit_dcc(r(1:150,:),dcc,max_iterations=25,tolerance=1.0e-5_dp)
   call assert_true(allocated(dcc%correlations).and.dcc%alpha>=0.0_dp.and.dcc%beta>=0.0_dp,'DCC fit output')
   call assert_true(dcc%alpha+dcc%beta<1.0_dp,'DCC stationarity')

   mean=[0.0_dp,0.0_dp]
   c=0.0_dp;c(1,1)=0.35_dp;c(2,1)=0.08_dp;c(2,2)=0.30_dp
   a=0.0_dp;a(1,1)=0.15_dp;a(2,2)=0.12_dp
   b=0.0_dp;b(1,1)=0.75_dp;b(2,2)=0.72_dp
   call bekk11_filter(r,mean,c,a,b,bcov,bres,istat)
   call assert_true(istat==mts_success.and.all(shape(bcov)==[k,k,n]),'BEKK filter')
   call assert_finite([bekk11_log_likelihood(r,mean,c,a,b)],'BEKK likelihood')

   call modified_cholesky_volatility(r,40,0.97_dp,mcov,coef,comp,istat)
   call assert_true(istat==mts_success.and.size(mcov,3)==n-40,'modified Cholesky volatility')
   call common_volatility_components(r,3,loadings,eigenvalues,transformed,status=istat)
   call assert_true(istat==mts_success.and.all(shape(loadings)==[k,k]),'common volatility components')
   call standardized_residuals(bres,bcov,stdres,istat)
   call volatility_diagnostics(bres,bcov,3,level,squared,istat)
   call assert_true(istat==mts_success.and.level%status==mts_success.and.squared%status==mts_success,'volatility diagnostics')
   call constrained_group_correlation(r,[1,1],span=100,unconstrained=ucor,constrained=gcor,status=istat)
   call assert_true(istat==mts_success.and.all(shape(gcor)==[k,k]),'constrained group correlation')

   print '(a)','test_volatility: PASS'
end program test_volatility
