program test_skewhyperbolic
  use skewhyperbolic
  implicit none
  integer, parameter :: n=30000
  real(dp) :: d0, expected, p, q, x(n), sm, sv, tm, tv, r(2), m2, m3
  type(skewhyp_fit_result) :: fit
  integer :: fails
  fails=0

  ! beta=0 is a scaled Student-t density.
  d0=dskewhyp(0.0_dp,0.0_dp,1.0_dp,0.0_dp,10.0_dp)
  expected=exp(log_gamma(5.5_dp)-0.5_dp*log(acos(-1.0_dp))-log_gamma(5.0_dp))
  call check(abs(d0-expected)<2.0e-10_dp,'beta=0 density',fails)

  p=pskewhyp(0.0_dp,0.0_dp,1.0_dp,0.0_dp,10.0_dp,tol=2.0e-8_dp)
  call check(abs(p-0.5_dp)<2.0e-5_dp,'symmetric cdf',fails)
  q=qskewhyp(0.8_dp,0.2_dp,1.3_dp,0.7_dp,12.0_dp,tol=3.0e-7_dp)
  p=pskewhyp(q,0.2_dp,1.3_dp,0.7_dp,12.0_dp,tol=3.0e-7_dp)
  call check(abs(p-0.8_dp)<2.0e-4_dp,'cdf quantile inversion',fails)

  call rskewhyp(x,0.5_dp,1.2_dp,0.4_dp,14.0_dp)
  sm=sum(x)/real(n,dp)
  sv=sum((x-sm)**2)/real(n-1,dp)
  tm=skewhyp_mean(0.5_dp,1.2_dp,0.4_dp,14.0_dp)
  tv=skewhyp_var(0.5_dp,1.2_dp,0.4_dp,14.0_dp)
  call check(abs(sm-tm)<0.035_dp,'rng mean',fails)
  call check(abs(sv-tv)<0.05_dp,'rng variance',fails)

  m2=skewhyp_moment(2,0.5_dp,1.2_dp,0.4_dp,14.0_dp,about=tm)
  m3=skewhyp_moment(3,0.5_dp,1.2_dp,0.4_dp,14.0_dp,about=tm)
  call check(abs(m2-tv)<2.0e-10_dp,'moment variance identity',fails)
  call check(abs(m3/(tv**1.5_dp)-skewhyp_skew(0.5_dp,1.2_dp,0.4_dp,14.0_dp))<2.0e-8_dp, &
             'moment skewness identity',fails)

  call skewhyp_calc_range(0.0_dp,1.0_dp,0.0_dp,10.0_dp,r,tol=1.0e-5_dp)
  call check(abs(dskewhyp(r(1),0.0_dp,1.0_dp,0.0_dp,10.0_dp)-1.0e-5_dp)<3.0e-7_dp,'range left',fails)
  call check(abs(dskewhyp(r(2),0.0_dp,1.0_dp,0.0_dp,10.0_dp)-1.0e-5_dp)<3.0e-7_dp,'range right',fails)

  call skewhyp_fit(x(1:800),fit,maxiter=220)
  call check(fit%param(2)>0.0_dp .and. fit%param(4)>0.0_dp,'fit valid parameters',fails)
  call check(abs(skewhyp_mean(fit%param(1),fit%param(2),fit%param(3),fit%param(4))-sm)<0.3_dp,'fit location',fails)

  if(fails==0)then
    print '(a)','test_skewhyperbolic: PASS'
  else
    print '(a,i0)','test_skewhyperbolic: FAIL ',fails
    error stop 1
  endif
contains
  subroutine check(ok,name,fails)
    logical,intent(in)::ok
    character(*),intent(in)::name
    integer,intent(inout)::fails
    if(.not.ok)then
      print '(a,a)','FAIL: ',name;fails=fails+1
    endif
  end subroutine
end program
