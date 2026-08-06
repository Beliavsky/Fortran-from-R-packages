program test_random
  use extra_distr
  use test_support
  implicit none
  integer, parameter :: n=12000
  real(dp), allocatable :: x(:),y(:)
  integer, allocatable :: ix(:),iy(:)

  call seed_rng(123456)
  x=rlaplace(20,0.5_dp,1.2_dp)
  call seed_rng(123456)
  y=rlaplace(20,0.5_dp,1.2_dp)
  call assert_true(all(x==y),'seeded continuous reproducibility')

  call seed_rng(98765)
  ix=rzip(20,2.0_dp,0.2_dp)
  call seed_rng(98765)
  iy=rzip(20,2.0_dp,0.2_dp)
  call assert_true(all(ix==iy),'seeded discrete reproducibility')

  call seed_rng(111)
  x=rlaplace(n,1.5_dp,0.8_dp)
  call assert_true(abs(sum(x)/real(n,dp)-1.5_dp)<0.04_dp,'Laplace random mean')

  x=rrayleigh(n,2.0_dp)
  call assert_true(all(x>=0.0_dp),'Rayleigh random support')
  call assert_true(abs(sum(x*x)/real(n,dp)-8.0_dp)<0.25_dp,'Rayleigh second moment')

  x=rkumar(n,2.0_dp,3.0_dp)
  call assert_true(all(x>=0.0_dp.and.x<=1.0_dp),'Kumaraswamy random support')
  x=rtnorm(n,0.0_dp,1.0_dp,-1.0_dp,2.0_dp)
  call assert_true(all(x>=-1.0_dp.and.x<=2.0_dp),'truncated normal random support')
  x=rinvgamma(n,4.0_dp,3.0_dp)
  call assert_true(all(x>0.0_dp),'inverse gamma random support')
  call assert_true(abs(sum(x)/real(n,dp)-1.0_dp)<0.06_dp,'inverse gamma random mean')

  ix=rbern(n,0.3_dp)
  call assert_true(abs(real(sum(ix),dp)/real(n,dp)-0.3_dp)<0.02_dp,'Bernoulli random mean')
  ix=rdunif(n,-3,4)
  call assert_true(all(ix>=-3.and.ix<=4),'discrete uniform random support')
  ix=rsign(n)
  call assert_true(all(abs(ix)==1),'Rademacher random support')
  ix=rtbinom(n,10,0.4_dp,2,7)
  call assert_true(all(ix>=2.and.ix<=7),'truncated binomial random support')
  ix=rtpois(n,4.0_dp,2,11)
  call assert_true(all(ix>=2.and.ix<=11),'truncated Poisson random support')

  x=rgev(1000,0.0_dp,1.0_dp,0.1_dp)
  call assert_true(all(x==x),'GEV random finite/defined')
  x=rwald(1000,1.0_dp,2.0_dp)
  call assert_true(all(x>0.0_dp),'Wald random support')
  x=rtlambda(1000,0.2_dp)
  call assert_true(all(x==x),'Tukey lambda random defined')

  call finish_tests('random generation')
end program test_random
