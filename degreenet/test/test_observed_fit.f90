! SPDX-License-Identifier: GPL-3.0-or-later
program test_observed_fit
  use degreenet_kinds, only : dp
  use degreenet_models, only : MODEL_DP, MODEL_YULE
  use degreenet_observed_fit, only : observed_fit_result, fit_grouped_model, fit_rounded_model
  implicit none
  type(observed_fit_result)::g,r
  call fit_grouped_model(MODEL_DP,[1,1,1,2,2,3,4,5,5,6,7],1,[2.5_dp],g,[1.01_dp],[10.0_dp],maxit=1200)
  if(.not.g%converged.or.g%theta(1)<=1.0_dp)error stop 1
  call fit_rounded_model(MODEL_YULE,[1,1,2,3,5,10,20,30,50],1,1000,[3.0_dp],r,[1.01_dp],[20.0_dp],maxit=1200)
  if(.not.r%converged.or.r%theta(1)<=1.0_dp)error stop 1
  print *, 'test_observed_fit: PASS'
end program
