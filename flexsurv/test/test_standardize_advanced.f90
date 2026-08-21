program test_standardize_advanced
  use flexsurv_kinds, only : dp
  use flexsurv_distributions, only : dist_exponential
  use flexsurv_fit, only : flexsurv_spec, initialize_spec
  use flexsurv_standardize_advanced, only : standsurv_acsurvival, standsurv_achazard, &
    standsurv_acrmst, standsurv_acquantile, standsurv_delta_acsurvival, standsurv_ci
  use relsurv_ratetable, only : ratetable_type, make_ratetable
  implicit none
  type(flexsurv_spec) :: spec
  type(ratetable_type) :: tab
  type(standsurv_ci) :: ci
  integer :: fails
  integer :: dims(1), factor(1), ncuts(1)
  real(dp) :: cuts(1,1), rate(1), rx(1,1), theta(1), cov(1,1)
  real(dp) :: t, s, h, r, q, truth

  fails=0
  call initialize_spec(spec,dist_exponential,1)
  theta(1)=log(0.2_dp)
  dims=1;factor=1;ncuts=0;cuts=0.0_dp;rate=0.1_dp
  tab=make_ratetable(dims,factor,cuts,ncuts,rate)
  rx=1.0_dp;t=2.0_dp

  s=standsurv_acsurvival(spec,theta,t,tab,rx)
  truth=exp(-0.3_dp*t)
  call check_close('acsurvival',s,truth,2.0e-10_dp,fails)

  h=standsurv_achazard(spec,theta,t,tab,rx)
  call check_close('achazard',h,0.3_dp,2.0e-10_dp,fails)

  r=standsurv_acrmst(spec,theta,t,tab,rx,nquad=48)
  truth=(1.0_dp-exp(-0.3_dp*t))/0.3_dp
  call check_close('acrmst',r,truth,2.0e-8_dp,fails)

  q=standsurv_acquantile(spec,theta,0.5_dp,tab,rx)
  call check_close('acquantile',q,log(2.0_dp)/0.3_dp,2.0e-8_dp,fails)

  cov(1,1)=0.01_dp
  ci=standsurv_delta_acsurvival(spec,theta,cov,t,tab,rx)
  call check_close('delta estimate',ci%estimate,exp(-0.6_dp),2.0e-10_dp,fails)
  if(.not.(ci%se>0.0_dp.and.ci%lower<ci%estimate.and.ci%upper>ci%estimate))then
    print *, 'FAIL delta interval',ci%se,ci%lower,ci%upper
    fails=fails+1
  end if

  if(fails==0)then
    print *, 'test_standardize_advanced: PASS'
  else
    print *, 'test_standardize_advanced: FAIL',fails
    error stop 1
  end if
contains
  subroutine check_close(name,a,b,tol,fails)
    character(len=*),intent(in)::name
    real(dp),intent(in)::a,b,tol
    integer,intent(inout)::fails
    if(abs(a-b)>tol*(1.0_dp+abs(b)))then
      print *, 'FAIL ',trim(name),a,b
      fails=fails+1
    end if
  end subroutine check_close
end program test_standardize_advanced
