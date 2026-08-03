! SPDX-License-Identifier: GPL-3.0-only
! Native Fortran translation of the computational core of rvinecopulib.
program test_bicop_core
  use rvinecopulib
  implicit none
  type(bicop_model) :: m
  real(dp) :: u,v,p,c_expected,cval
  m=make_bicop(bicop_indep)
  call assert_close(m%pdf(0.2_dp,0.7_dp),1.0_dp,1.0e-14_dp,'independence pdf')
  call assert_close(m%cdf(0.2_dp,0.7_dp),0.14_dp,1.0e-14_dp,'independence cdf')
  m=make_bicop(bicop_clayton,0,[2.0_dp])
  c_expected=(0.5_dp**(-2.0_dp)+0.5_dp**(-2.0_dp)-1.0_dp)**(-0.5_dp)
  call assert_close(m%cdf(0.5_dp,0.5_dp),c_expected,2.0e-13_dp,'Clayton cdf')
  call assert_close(m%tau(),0.5_dp,1.0e-13_dp,'Clayton tau')
  u=0.3_dp; v=0.7_dp; p=m%hfunc1(u,v)
  call assert_close(m%hinv1(u,p),v,2.0e-12_dp,'h inverse 1')
  p=m%hfunc2(u,v)
  call assert_close(m%hinv2(p,v),u,2.0e-12_dp,'h inverse 2')
  m=make_bicop(bicop_clayton,90,[2.0_dp])
  if (m%tau()>=0.0_dp) error stop 'rotation tau sign'
  if (m%pdf(u,v)<=0.0_dp) error stop 'rotation pdf'
  m=make_bicop(bicop_gaussian,0,[0.5_dp])
  call assert_close(m%tau(),2.0_dp/pi*asin(0.5_dp),1.0e-13_dp,'Gaussian tau')
  cval=m%cdf(u,v)
  if (cval<=max(0.0_dp,u+v-1.0_dp) .or. cval>=min(u,v)) error stop 'Gaussian cdf bounds'
  print '(a)', 'test_bicop_core: PASS'
contains
  subroutine assert_close(a,b,tol,msg)
    real(dp),intent(in)::a,b,tol
    character(len=*),intent(in)::msg
    if(abs(a-b)>tol) then
      print *,msg,a,b
      error stop 1
    end if
  end subroutine
end program
