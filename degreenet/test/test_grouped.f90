! SPDX-License-Identifier: GPL-3.0-or-later
program test_grouped
  use degreenet_kinds, only : dp
  use degreenet_models, only : MODEL_DP, grouped_probability, model_pmf
  use degreenet_observation, only : rounded_probability, rounded_bin, grouped_loglik
  use degreenet_diagnostics, only : mands_result, modified_anderson_darling
  implicit none
  integer::lo,hi,k
  real(dp)::p,s
  type(mands_result)::mr
  call rounded_bin(20,lo,hi)
  if(lo/=16.or.hi/=24)error stop 1
  p=rounded_probability(MODEL_DP,[3.0_dp],20,1);s=0.0_dp
  do k=16,24;s=s+model_pmf(MODEL_DP,[3.0_dp],k,1);end do
  if(abs(p-s)>1e-14_dp)error stop 1
  p=grouped_probability(MODEL_DP,[3.0_dp],5,1);s=0.0_dp
  do k=5,10;s=s+model_pmf(MODEL_DP,[3.0_dp],k,1);end do
  if(abs(p-s)>1e-14_dp)error stop 1
  if(grouped_loglik(MODEL_DP,[3.0_dp],[1,2,5,6],1)>=0.0_dp)error stop 1
  call modified_anderson_darling(MODEL_DP,[3.0_dp],[1,1,1,2,2,3,4,5],1,1000,mr)
  if(mr%statistic<0.0_dp)error stop 1
  print *, 'test_grouped: PASS'
end program
