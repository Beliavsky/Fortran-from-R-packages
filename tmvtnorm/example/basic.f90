program basic
  use tmvtnorm
  implicit none
  real(dp)::mu(2),sigma(2,2),lower(2),upper(2)
  real(dp),allocatable::x(:,:)
  type(tmvnorm_moments_t)::mom
  mu=[0.0_dp,0.0_dp]
  sigma=reshape([1.0_dp,0.6_dp,0.6_dp,1.5_dp],[2,2])
  lower=[-1.0_dp,-0.5_dp]
  upper=[1.5_dp,2.0_dp]
  call mtmvnorm(mu,sigma,lower,upper,mom)
  print '(a,2f12.6)','truncated mean: ',mom%mean
  x=rtmvnorm(5,mu,sigma,lower,upper,algorithm=algorithm_gibbs,burnin=20,seed=1234)
  print '(a)','five Gibbs draws:'
  print '(2f12.6)',transpose(x)
end program
