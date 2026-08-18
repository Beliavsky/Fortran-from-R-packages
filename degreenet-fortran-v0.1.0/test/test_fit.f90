! SPDX-License-Identifier: GPL-3.0-or-later
program test_fit
  use degreenet_kinds, only : dp
  use degreenet_models, only : MODEL_GEOM, MODEL_DP, loglik_model
  use degreenet_fit, only : fit_result, fit_degree_model
  implicit none
  integer, parameter :: x(20)=[1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,3,3,4,5,6]
  type(fit_result)::fit
  real(dp)::ll1,ll2
  call fit_degree_model(MODEL_GEOM,x,1,1000,[2.0_dp],fit,[1.001_dp],[100.0_dp],maxit=1000)
  if(.not.fit%converged)then;print *,'FAIL fit convergence';error stop 1;end if
  if(fit%theta(1)<=1.0_dp.or.fit%theta(1)>5.0_dp)then;print *,'FAIL fit theta',fit%theta;error stop 1;end if
  ll1=loglik_model(MODEL_GEOM,[2.0_dp],x,1,1000)
  ll2=loglik_model(MODEL_GEOM,fit%theta,x,1,1000)
  if(ll2<ll1-1e-8_dp)then;print *,'FAIL fit likelihood';error stop 1;end if
  if(loglik_model(MODEL_DP,[3.0_dp],x,1,1000)>0.0_dp)error stop 1
  print *, 'test_fit: PASS'
end program
