! SPDX-License-Identifier: Artistic-2.0
program factor_vecm_analysis
   use mts
   implicit none
   integer,parameter::n=400,k=4
   real(dp)::panel(n,k),common,eps(k),cointegrated(n,2),beta(2,1)
   type(factor_result)::pca
   type(vecm_model)::vecm
   integer::t

   call set_random_seed(4404)
   cointegrated(1,:)=[0.0_dp,0.0_dp]
   call random_multivariate_normal([0.0_dp,0.0_dp,0.0_dp,0.0_dp],eye(k),eps)
   common=sin(0.03_dp)+0.2_dp*eps(1)
   panel(1,:)=[1.0_dp,0.8_dp,-0.5_dp,0.3_dp]*common+0.15_dp*eps
   do t=2,n
      call random_multivariate_normal([0.0_dp,0.0_dp,0.0_dp,0.0_dp],eye(k),eps)
      common=sin(0.03_dp*real(t,dp))+0.2_dp*eps(1)
      panel(t,:)=[1.0_dp,0.8_dp,-0.5_dp,0.3_dp]*common+0.15_dp*eps
      cointegrated(t,1)=cointegrated(t-1,1)+0.1_dp*eps(1)
      cointegrated(t,2)=cointegrated(t,1)+0.5_dp*(cointegrated(t-1,2)-cointegrated(t-1,1))+0.05_dp*eps(2)
   end do
   call principal_components(panel,1,pca,standardize=.true.)
   beta(:,1)=[1.0_dp,-1.0_dp]
   call fit_vecm_known_beta(cointegrated,1,beta,vecm,include_constant=.true.)

   print '(a,4(1x,f9.5))','first factor loadings:',pca%loadings(:,1)
   print '(a,2(1x,f9.5))','VECM adjustment coefficients:',vecm%alpha(:,1)
end program factor_vecm_analysis
