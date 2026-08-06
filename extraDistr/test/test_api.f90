program test_api
  use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
  use extra_distr
  use test_support
  implicit none
  real(dp) :: p,d
  integer :: q
  integer, allocatable :: z(:)
  real(dp), allocatable :: rz(:)

  d=dlaplace(0.3_dp,0.0_dp,1.2_dp)
  call assert_close(exp(dlaplace(0.3_dp,0.0_dp,1.2_dp,log_p=.true.)),d,1.0e-13_dp,'continuous log density')
  p=plaplace(0.3_dp,0.0_dp,1.2_dp)
  call assert_close(plaplace(0.3_dp,0.0_dp,1.2_dp,lower_tail=.false.),1.0_dp-p,1.0e-13_dp,'continuous upper tail')
  call assert_close(exp(plaplace(0.3_dp,0.0_dp,1.2_dp,log_p=.true.)),p,1.0e-13_dp,'continuous log CDF')
  call assert_close(qlaplace(log(p),0.0_dp,1.2_dp,log_p=.true.),0.3_dp,1.0e-12_dp,'log probability quantile')
  call assert_close(qlaplace(1.0_dp-p,0.0_dp,1.2_dp,lower_tail=.false.),0.3_dp,1.0e-12_dp,'upper-tail quantile')

  q=qzip(0.72_dp,2.3_dp,0.2_dp)
  call assert_true(pzip(q,2.3_dp,0.2_dp)>=0.72_dp.and.(q==0.or.pzip(q-1,2.3_dp,0.2_dp)<0.72_dp),'ZIP quantile convention')
  q=qtbinom(0.6_dp,9,0.4_dp,2,7)
  call assert_true(ptbinom(q,9,0.4_dp,2,7)>=0.6_dp.and.(q==2.or.ptbinom(q-1,9,0.4_dp,2,7)<0.6_dp),'truncated-binomial quantile convention')

  call assert_true(ieee_is_nan(dbetapr(1.0_dp,-1.0_dp,2.0_dp)),'invalid beta-prime parameters')
  call assert_true(ieee_is_nan(ddunif(0,2,1)),'invalid discrete-uniform bounds')
  call assert_true(ieee_is_nan(dbvnorm(0.0_dp,0.0_dp,cor=1.1_dp)),'invalid correlation')
  call assert_close(ddgamma(2,3.0_dp,rate=2.0_dp),ddgamma(2,3.0_dp,scale=0.5_dp),1.0e-14_dp,'discrete gamma rate alias')
  call assert_close(dgpois(2,3.0_dp,rate=2.0_dp),dgpois(2,3.0_dp,scale=0.5_dp),1.0e-14_dp,'gamma-Poisson rate alias')
  call assert_true(ieee_is_nan(pbern(0,-0.2_dp)),'invalid Bernoulli CDF parameter')
  call assert_true(ieee_is_nan(pbbinom(1,-1,2.0_dp,3.0_dp)),'invalid beta-binomial CDF parameter')
  call assert_true(ieee_is_nan(pbnbinom(1,-1.0_dp,2.0_dp,3.0_dp)),'invalid beta-negative-binomial CDF parameter')
  call assert_true(ieee_is_nan(pcat(1,[-1.0_dp,2.0_dp])),'invalid categorical CDF parameter')
  call assert_true(ieee_is_nan(pdgamma(1,-2.0_dp)),'invalid discrete-gamma CDF parameter')
  call assert_true(ieee_is_nan(pdnorm(1,sd=-1.0_dp)),'invalid discrete-normal CDF parameter')
  call assert_true(ieee_is_nan(pdunif(1,3,2)),'invalid discrete-uniform CDF bounds')
  call assert_true(ieee_is_nan(pdweibull(1,1.2_dp,2.0_dp)),'invalid discrete-Weibull CDF parameter')
  call assert_true(ieee_is_nan(pgpois(1,-1.0_dp,rate=1.0_dp)),'invalid gamma-Poisson CDF parameter')
  call assert_true(ieee_is_nan(plgser(1,1.2_dp)),'invalid log-series CDF parameter')
  call assert_true(ieee_is_nan(pnhyper(2,3,2,3)),'invalid negative-hypergeometric CDF parameter')
  call assert_true(ieee_is_nan(dtbinom(2,5,-0.1_dp,1,4)),'invalid truncated-binomial parameter')
  call assert_true(ieee_is_nan(dtpois(2,-1.0_dp,1,4)),'invalid truncated-Poisson parameter')
  call assert_true(ieee_is_nan(dzip(1,-1.0_dp,0.2_dp)),'invalid ZIP parameter')
  call assert_true(ieee_is_nan(dzib(1,5,0.4_dp,1.2_dp)),'invalid ZIB parameter')
  call assert_true(ieee_is_nan(dzinb(1,3.0_dp,-0.4_dp,0.2_dp)),'invalid ZINB parameter')
  call assert_true(ieee_is_nan(ddirmnom([1,1],2,[-1.0_dp,2.0_dp])),'invalid Dirichlet-multinomial parameter')
  call assert_true(ieee_is_nan(dmnom([1,1],2,[-1.0_dp,2.0_dp])),'invalid multinomial parameter')
  call assert_true(ieee_is_nan(dmvhyper([1,1],[-1,3],2)),'invalid multivariate-hypergeometric parameter')
  call assert_true(ieee_is_nan(dmixnorm(0.0_dp,[0.0_dp],[1.0_dp],[-1.0_dp])),'invalid normal-mixture parameter')
  call assert_true(ieee_is_nan(dmixpois(0,[1.0_dp],[-1.0_dp])),'invalid Poisson-mixture parameter')


  call assert_close(dslash(2.0_dp,2.0_dp,2.0_dp),1.0_dp/(2.0_dp*sqrt(2.0_dp*acos(-1.0_dp))),1.0e-13_dp,'source slash center')
  call assert_close(dslash(2.0_dp,2.0_dp,2.0_dp,source_compatible=.false.),1.0_dp/(4.0_dp*sqrt(2.0_dp*acos(-1.0_dp))),1.0e-13_dp,'corrected slash center')
  call assert_true(pdlaplace(0,2.0_dp,0.5_dp)/=pdlaplace(0,2.0_dp,0.5_dp,source_compatible=.false.),'discrete Laplace compatibility switch')
  call assert_true(ddirichlet([0.2_dp,0.2_dp],[1.0_dp,1.0_dp])>0.0_dp,'source Dirichlet off-simplex behavior')
  call assert_close(ddirichlet([0.2_dp,0.2_dp],[1.0_dp,1.0_dp],source_compatible=.false.),0.0_dp,0.0_dp,'corrected Dirichlet support')
  call seed_rng(19)
  rz=rdlaplace(20,0.5_dp,0.4_dp)
  call assert_true(all(abs(rz-0.5_dp-real(nint(rz-0.5_dp),dp))<1.0e-12_dp),'shifted discrete-Laplace random support')

  z=rbern(0,0.3_dp)
  call assert_int(size(z),0,'zero-length random result')
  call assert_close(dpareto(0.5_dp,2.0_dp,1.0_dp),0.0_dp,0.0_dp,'Pareto support')
  call assert_close(dtnorm(-2.0_dp,0.0_dp,1.0_dp,-1.0_dp,1.0_dp),0.0_dp,0.0_dp,'truncated normal support')
  call assert_close(ddweibull(-1,0.3_dp,1.2_dp),0.0_dp,0.0_dp,'discrete Weibull support')

  call finish_tests('API tails, logs, and errors')
end program test_api
