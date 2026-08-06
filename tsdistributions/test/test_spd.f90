program test_spd
  use tsdistributions
  implicit none
  type(rng_state)::rng
  type(distribution_parameters)::p
  type(spd_specification)::spec
  type(spd_fit)::fit
  real(dp),allocatable::x(:),z(:)
  real(dp)::pr,q
  integer::i
  call seed_rng(rng,77123_i8)
  p=distribution_parameters(mu=0.0_dp,sigma=1.0_dp,shape=6.0_dp)
  x=rdist('std',800,rng,p)
  spec=spd_modelspec(x,0.1_dp,0.9_dp,'normal')
  fit=estimate_spd(spec)
  call assert_true(fit%status==tsd_success,'SPD fit')
  do i=1,9
    pr=real(i,dp)/10.0_dp;q=qspd(pr,fit)
    call assert_true(abs(pspd(q,fit)-pr)<2.0e-5_dp,'SPD CDF/quantile')
    call assert_true(dspd(q,fit)>0.0_dp,'SPD density')
  end do
  z=rspd(100,rng=rng,fit=fit)
  call assert_true(size(z)==100.and.all(abs(z)<huge(1.0_dp)),'SPD random generation')
  print '(a)','test_spd: PASS'
contains
  subroutine assert_true(condition,message)
    logical,intent(in)::condition;character(len=*),intent(in)::message
    if(.not.condition)then;write(*,'(a)')'FAIL: '//message;error stop 1;end if
  end subroutine
end program
