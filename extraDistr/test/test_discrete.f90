program test_discrete
  use extra_distr
  use test_support
  implicit none
  integer :: k
  real(dp) :: s, p
  real(dp), parameter :: tol = 2.0e-8_dp
  real(dp) :: probs(3)

  call assert_close(dbern(0,0.3_dp)+dbern(1,0.3_dp),1.0_dp,tol,'Bernoulli normalization')
  call assert_int(qbern(0.69_dp,0.3_dp),0,'Bernoulli quantile low')
  call assert_int(qbern(0.71_dp,0.3_dp),1,'Bernoulli quantile high')

  s=0.0_dp; do k=0,8; s=s+dbbinom(k,8,2.0_dp,3.0_dp); end do
  call assert_close(s,1.0_dp,tol,'beta-binomial normalization')
  call assert_close(pbbinom(3,8,2.0_dp,3.0_dp),sum([(dbbinom(k,8,2.0_dp,3.0_dp),k=0,3)]),tol,'beta-binomial CDF')

  s=0.0_dp; do k=0,1000; s=s+dbnbinom(k,2.5_dp,4.0_dp,5.0_dp); end do
  call assert_close(s,1.0_dp,2.0e-6_dp,'beta-negative-binomial normalization')

  probs=[1.0_dp,2.0_dp,3.0_dp]
  call assert_close(sum([(dcat(k,probs),k=1,3)]),1.0_dp,tol,'categorical normalization')
  call assert_int(qcat(0.51_dp,probs),3,'categorical quantile')

  s=0.0_dp; do k=0,100; s=s+ddgamma(k,2.0_dp,1.5_dp); end do
  call assert_close(s,1.0_dp,1.0e-8_dp,'discrete gamma normalization')
  s=0.0_dp; do k=-100,100; s=s+ddlaplace(k,0.0_dp,0.4_dp); end do
  call assert_close(s,1.0_dp,1.0e-10_dp,'discrete Laplace normalization')
  s=0.0_dp; do k=-30,30; s=s+ddnorm(k,0.3_dp,1.7_dp); end do
  call assert_close(s,1.0_dp,1.0e-10_dp,'discrete normal normalization')
  call assert_close(sum([(ddunif(k,-2,4),k=-2,4)]),1.0_dp,tol,'discrete uniform normalization')
  call assert_int(qdunif(0.5_dp,-2,4),1,'discrete uniform median')

  s=0.0_dp; do k=0,200; s=s+ddweibull(k,0.2_dp,1.4_dp); end do
  call assert_close(s,1.0_dp,1.0e-9_dp,'discrete Weibull normalization')
  p=0.63_dp; k=qdweibull(p,0.2_dp,1.4_dp)
  call assert_true(pdweibull(k,0.2_dp,1.4_dp)>=p.and.(k==0.or.pdweibull(k-1,0.2_dp,1.4_dp)<p),'discrete Weibull quantile')

  s=0.0_dp; do k=0,250; s=s+dgpois(k,2.3_dp,1.4_dp); end do
  call assert_close(s,1.0_dp,1.0e-9_dp,'gamma-Poisson normalization')
  s=0.0_dp; do k=1,5000; s=s+dlgser(k,0.6_dp); end do
  call assert_close(s,1.0_dp,1.0e-9_dp,'log-series normalization')
  p=0.7_dp; k=qlgser(p,0.6_dp)
  call assert_true(plgser(k,0.6_dp)>=p.and.(k==1.or.plgser(k-1,0.6_dp)<p),'log-series quantile')

  s=0.0_dp; do k=3,13; s=s+dnhyper(k,10,8,3); end do
  call assert_close(s,1.0_dp,1.0e-10_dp,'negative hypergeometric normalization')
  p=0.55_dp; k=qnhyper(p,10,8,3)
  call assert_true(pnhyper(k,10,8,3)>=p.and.(k==3.or.pnhyper(k-1,10,8,3)<p),'negative hypergeometric quantile')

  s=0.0_dp; do k=-40,40; s=s+dskellam(k,2.0_dp,3.0_dp); end do
  call assert_close(s,1.0_dp,2.0e-8_dp,'Skellam normalization')

  s=0.0_dp; do k=2,7; s=s+dtbinom(k,10,0.35_dp,2,7); end do
  call assert_close(s,1.0_dp,tol,'truncated binomial normalization')
  s=0.0_dp; do k=2,12; s=s+dtpois(k,4.0_dp,2,12); end do
  call assert_close(s,1.0_dp,tol,'truncated Poisson normalization')

  s=0.0_dp; do k=0,10; s=s+dzib(k,10,0.3_dp,0.25_dp); end do
  call assert_close(s,1.0_dp,tol,'zero-inflated binomial normalization')
  s=0.0_dp; do k=0,500; s=s+dzinb(k,3.0_dp,0.55_dp,0.2_dp); end do
  call assert_close(s,1.0_dp,1.0e-9_dp,'zero-inflated negative-binomial normalization')
  s=0.0_dp; do k=0,100; s=s+dzip(k,3.5_dp,0.3_dp); end do
  call assert_close(s,1.0_dp,1.0e-10_dp,'zero-inflated Poisson normalization')

  call assert_close(exp(dzip(2,3.5_dp,0.3_dp,log_p=.true.)),dzip(2,3.5_dp,0.3_dp),tol,'discrete log density')
  call assert_close(pzip(2,3.5_dp,0.3_dp,lower_tail=.false.),1.0_dp-pzip(2,3.5_dp,0.3_dp),tol,'discrete upper tail')

  call finish_tests('discrete distributions')
end program test_discrete
