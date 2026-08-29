program test_estimation
  use tmvtnorm
  implicit none
  real(dp)::mu(1),s(1,1),lo(1),up(1),mu0(1),s0(1,1)
  real(dp),allocatable::x(:,:),g(:,:)
  type(tmvnorm_fit_t)::fit
  type(tmvnorm_gmm_fit_t)::gf
  mu=[0.4_dp]
  s(1,1)=1.1_dp
  lo=[-1.0_dp]
  up=[2.0_dp]
  x=rtmvnorm_gibbs(300,mu,s,lo,up,seed=2026)
  mu0=[0.0_dp]
  s0(1,1)=1.0_dp
  call mle_tmvnorm(x,lo,up,mu0,s0,fit,cholesky=.true.,maxit=1200,tol=1e-6_dp)
  if(abs(fit%mean(1)-mu(1))>0.35_dp) then
  print *,'mle mean',fit%mean
  error stop 1
  end if
  if(abs(fit%sigma(1,1)-s(1,1))>0.45_dp) then
  print *,'mle var',fit%sigma
  error stop 1
  end if
  g=gmm_moments_manjunath_wilhelm(x,mu,s,lo,up)
  if(maxval(abs(sum(g,dim=1)/real(size(g,1),dp)))>0.15_dp) error stop 'MW moments too far from zero'
  call gmm_tmvnorm(x,lo,up,mu0,s0,gf,method=gmm_mw,cholesky=.true.,maxit=800,tol=2e-5_dp)
  if(abs(gf%mean(1)-mu(1))>0.45_dp) then
  print *,'gmm mean',gf%mean
  error stop 1
  end if
  print *, 'test_estimation: ok'
end program
