! SPDX-License-Identifier: GPL-2.0-or-later
! This file is part of a modern Fortran translation of PerformanceAnalytics.
! Copyright (C) original PerformanceAnalytics authors and translation contributors.
! Distributed under GNU GPL version 2 or, at your option, any later version.
program test_capm_moments
  use kinds_mod, only: dp
  use capm_mod
  use comoments_mod
  use portfolio_risk_mod
  use statistics_mod, only: mean_value, variance_value, centered_moment
  use test_support_mod
  implicit none
  real(dp) :: market(8), asset(8), rf(8), r(6,2), cov(2,2), m3(2,4), m4(2,8)
  real(dp) :: v3(4), v3b(4), v4(5), v4b(5), w(2), mu(2), pr(6)
  real(dp) :: comp(2), total, ewma(2,2)
  type(sfm_result) :: fit
  type(market_timing_result) :: timing
  integer :: i

  market=[-0.04_dp,-0.02_dp,-0.01_dp,0.0_dp,0.01_dp,0.02_dp,0.03_dp,0.05_dp]
  asset=0.001_dp+1.5_dp*market
  rf=0.0_dp
  call sfm_fit(asset,market,rf,fit)
  call assert_close(fit%alpha,0.001_dp,1.0e-12_dp,'sfm alpha')
  call assert_close(fit%beta,1.5_dp,1.0e-12_dp,'sfm beta')
  call assert_close(fit%r_squared,1.0_dp,1.0e-12_dp,'sfm r2')
  call assert_close(capm_beta_bull(asset,market,rf),1.5_dp,1.0e-12_dp,'bull beta')
  call assert_close(capm_beta_bear(asset,market,rf),1.5_dp,1.0e-12_dp,'bear beta')
  call assert_true(specific_risk(asset,market,rf)<1.0e-12_dp,'specific risk')

  asset=0.002_dp+1.1_dp*market+0.8_dp*market*market
  call treynor_mazuy_fit(asset,market,rf,timing)
  call assert_close(timing%gamma,0.8_dp,1.0e-10_dp,'timing gamma')

  r(:,1)=[-0.02_dp,-0.01_dp,0.00_dp,0.01_dp,0.02_dp,0.03_dp]
  r(:,2)=[0.03_dp,0.01_dp,-0.01_dp,0.00_dp,0.02_dp,0.04_dp]
  call covariance_matrix(r,cov,.false.)
  call assert_close(cov(1,1),variance_value(r(:,1),.false.),1.0e-13_dp,'cov variance')
  call coskewness_matrix(r,m3)
  call cokurtosis_matrix(r,m4)
  call coskewness_unique(r,v3)
  call cokurtosis_unique(r,v4)
  call m3_mat_to_vec(m3,2,v3b)
  call m4_mat_to_vec(m4,2,v4b)
  call assert_vector_close(v3,v3b,1.0e-13_dp,'m3 conversion')
  call assert_vector_close(v4,v4b,1.0e-13_dp,'m4 conversion')
  call assert_true(m3_inner_product(v3,v3,2)>=0.0_dp,'m3 norm')
  call assert_true(m4_inner_product(v4,v4,2)>=0.0_dp,'m4 norm')

  w=[0.6_dp,0.4_dp];mu=[mean_value(r(:,1)),mean_value(r(:,2))]
  pr=matmul(r,w)
  call assert_close(portfolio_mean(w,mu),mean_value(pr),1.0e-13_dp,'portfolio mean')
  call assert_close(portfolio_variance(w,cov),variance_value(pr,.false.),1.0e-13_dp,'portfolio variance')
  call assert_close(portfolio_third_moment(w,m3),centered_moment(pr,3),1.0e-13_dp,'portfolio third')
  call assert_close(portfolio_fourth_moment(w,m4),centered_moment(pr,4),1.0e-13_dp,'portfolio fourth')
  call assert_true(portfolio_gaussian_es(w,mu,cov,0.95_dp)>portfolio_gaussian_var(w,mu,cov,0.95_dp),'portfolio es')
  call portfolio_risk_contributions(w,mu,cov,m3,m4,0.95_dp,'modified_var',comp,total)
  call assert_true(all(abs(comp)<100.0_dp),'portfolio risk contributions')
  call assert_close(herfindahl_index([0.5_dp,0.5_dp]),0.5_dp,1.0e-13_dp,'herfindahl')
  call ewma_covariance(r,0.94_dp,ewma)
  call assert_true(all(abs(ewma-transpose(ewma))<1.0e-12_dp),'ewma symmetry')
  do i=1,2
    call assert_true(ewma(i,i)>=0.0_dp,'ewma diagonal')
  end do
  write(*,'(a)')'CAPM, co-moment, and portfolio-risk tests passed.'
end program test_capm_moments
