program test_core
  use tolerance
  implicit none
  real(dp) :: v
  integer :: fail,i
  fail=0
  call chk(normal_cdf(0.0_dp),0.5_dp,2e-14_dp,'normal cdf')
  call chk(normal_quantile(0.975_dp),1.959963984540054_dp,2e-10_dp,'normal quantile')
  call chk(polygamma(1,1.0_dp),1.6449340668482264_dp,2e-10_dp,'trigamma')
  call chk(beta_i(0.5_dp,2.0_dp,3.0_dp),0.6875_dp,2e-12_dp,'beta cdf')
  v=sum([(dpoislind(i,2.0_dp),i=0,100)])
  call chk(v,1.0_dp,1e-12_dp,'poisson-lindley normalization')
  if(qnhyper(0.5_dp,20,30,5)<0)call bad('negative hypergeometric quantile')
  if(fail==0)then;print '(a)','test_core: PASS';else;error stop 1;end if
contains
  subroutine chk(x,y,tol,name)
    real(dp),intent(in)::x,y,tol;character(len=*),intent(in)::name
    if(abs(x-y)>tol)then;print *,trim(name),x,y;fail=fail+1;end if
  end subroutine
  subroutine bad(name)
    character(len=*),intent(in)::name;print *,trim(name);fail=fail+1
  end subroutine
end program
