program test_random_moments
  use tsdistributions
  implicit none
  type(rng_state)::rng
  type(distribution_parameters)::p
  type(moment_summary)::m
  real(dp),allocatable::x(:)
  real(dp)::out(4),sk0,ku0
  integer::status
  call seed_rng(rng,20260804_i8)
  p=distribution_parameters(mu=1.0_dp,sigma=2.0_dp,skew=1.0_dp,shape=8.0_dp)
  x=rdist('std',20000,rng,p)
  call assert_true(abs(sum(x)/real(size(x),dp)-1.0_dp)<0.08_dp,'Student random mean')
  call assert_true(abs(sample_sd(x)-2.0_dp)<0.10_dp,'Student random SD')
  m=distribution_moments('std',shape=8.0_dp)
  call assert_true(abs(m%skewness)<1.0e-14_dp.and.abs(m%excess_kurtosis-1.5_dp)<1.0e-12_dp,'Student moments')
  sk0=dskewness('norm');ku0=dkurtosis('norm')
  call assert_true(abs(sk0)<1.0e-14_dp.and.abs(ku0)<1.0e-14_dp,'normal moments')
  call nigtransform(0.0_dp,1.0_dp,0.2_dp,1.5_dp,out,status)
  call assert_true(status==tsd_success.and.all(abs(out)<huge(1.0_dp)),'NIG transform')
  print '(a)','test_random_moments: PASS'
contains
  subroutine assert_true(condition,message)
    logical,intent(in)::condition;character(len=*),intent(in)::message
    if(.not.condition)then;write(*,'(a)')'FAIL: '//message;error stop 1;end if
  end subroutine
end program
