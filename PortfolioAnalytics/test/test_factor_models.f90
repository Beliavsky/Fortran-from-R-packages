! SPDX-License-Identifier: GPL-3.0-only
program test_factor_models
  use portfolio_analytics
  use test_support
  implicit none
  real(dp) :: r(8,3),sample_mu(3),sample_sigma(3,3),model_sigma(3,3)
  real(dp) :: beta(2),stock2(2),stock3(2),stock4(2),sigma2(2,2),m3(2,4),m4(2,8)
  real(dp) :: clean(8,3),shr_mu(3),shr_sigma(3,3),rob_mu(3),rob_sigma(3,3)
  type(factor_model_result) :: model
  integer :: i,info

  do i=1,8
    r(i,1)=0.01_dp*real(i-4,dp)
    r(i,2)=2.0_dp*r(i,1)+0.001_dp*(-1.0_dp)**i
    r(i,3)=-r(i,1)+0.002_dp*(-1.0_dp)**(i+1)
  end do
  call sample_moments(r,sample_mu,sample_sigma)
  call fit_statistical_factor_model(r,3,model,info)
  call assert_true(info==0,'factor model fit')
  call assert_true(model%eigenvalues(1)>=model%eigenvalues(2),'ordered PCA eigenvalues')
  call factor_model_covariance(model,model_sigma)
  call assert_true(maxval(abs(model_sigma-sample_sigma))<1.0e-9_dp,'full PCA covariance reconstruction')

  beta=[1.0_dp,2.0_dp]
  stock2=[0.1_dp,0.2_dp]
  stock3=[0.01_dp,0.02_dp]
  stock4=[0.04_dp,0.08_dp]
  call covariance_sf(beta,stock2,0.5_dp,sigma2)
  call assert_close(sigma2(1,2),1.0_dp,1.0e-14_dp,'single-factor covariance')
  call coskewness_sf(beta,stock3,0.25_dp,m3)
  call assert_close(m3(1,1),0.26_dp,1.0e-14_dp,'single-factor coskewness diagonal')
  call cokurtosis_sf([0.0_dp,0.0_dp],stock2,stock4,0.5_dp,0.75_dp,m4)
  call assert_close(m4(1,1),stock4(1),1.0e-14_dp,'residual fourth moment')

  call winsorize_returns(r,0.1_dp,0.9_dp,clean)
  call covariance_shrinkage(r,1.0_dp,shr_mu,shr_sigma)
  call assert_close(shr_sigma(1,2),0.0_dp,1.0e-14_dp,'diagonal shrinkage target')
  call robust_covariance_huber(r,2.5_dp,50,rob_mu,rob_sigma)
  call assert_true(all([(rob_sigma(i,i)>0.0_dp,i=1,3)]),'robust covariance positive diagonal')
  print '(a)','test_factor_models: PASS'
end program test_factor_models
