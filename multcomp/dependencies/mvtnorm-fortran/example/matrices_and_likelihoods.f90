! SPDX-License-Identifier: GPL-2.0-only
program matrices_and_likelihoods
  use mvtnorm
  implicit none
  real(dp) :: mean(3),sigma(3,3),observations(2,3),given(1)
  integer :: which(1)
  real(dp),allocatable :: chol(:,:),invchol(:,:),correlation(:,:)
  type(conditional_result) :: conditional
  type(likelihood_result) :: likelihood
  logical :: ok
  character(len=256) :: message

  mean=[0.1_dp,-0.2_dp,0.3_dp]
  sigma=reshape([1.0_dp,0.3_dp,0.1_dp, 0.3_dp,1.2_dp,-0.2_dp, 0.1_dp,-0.2_dp,0.9_dp],[3,3])
  observations=reshape([0.2_dp,-0.4_dp, -0.1_dp,0.5_dp, 0.7_dp,0.0_dp],[2,3])

  chol=cov2chol(sigma,ok,message)
  if(.not.ok) error stop trim(message)
  invchol=chol2invchol(chol,ok)
  correlation=invchol2cor(invchol)
  likelihood=ldmvnorm(observations,mean,sigma,.true.)
  which=[1]; given=[0.5_dp]
  conditional=conditional_mvnormal(mean,sigma,which,given)

  print '(a,f12.6)', 'total log likelihood: ',likelihood%total
  print '(a,3f10.5)', 'correlation row 1: ',correlation(1,:)
  print '(a,2f10.5)', 'conditional mean: ',conditional%mean
end program matrices_and_likelihoods
