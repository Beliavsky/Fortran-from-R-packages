! SPDX-License-Identifier: GPL-3.0-only
! Native Fortran translation of the computational core of rvinecopulib.
program test_cvine
  use rvinecopulib
  implicit none
  type(cvine_model) :: vine,indep
  real(dp) :: x(4,1),z(4,1),back(4,1)
  vine=make_cvine(4)
  vine%pair(1,2)=make_bicop(bicop_gaussian,0,[0.4_dp])
  vine%pair(1,3)=make_bicop(bicop_clayton,0,[1.0_dp])
  vine%pair(1,4)=make_bicop(bicop_gumbel,0,[1.3_dp])
  vine%pair(2,3)=make_bicop(bicop_frank,0,[2.0_dp])
  vine%pair(2,4)=make_bicop(bicop_gaussian,0,[-0.2_dp])
  vine%pair(3,4)=make_bicop(bicop_joe,0,[1.2_dp])
  x(:,1)=[0.2_dp,0.45_dp,0.7_dp,0.8_dp]
  call vine%rosenblatt(x,z)
  call vine%inverse_rosenblatt(z,back)
  if (maxval(abs(back-x))>2.0e-10_dp) error stop 'C-vine Rosenblatt roundtrip'
  indep=make_cvine(4)
  call assert_close(indep%pdf(x(:,1)),1.0_dp,1.0e-13_dp,'independent C-vine')
  call vine%truncate(1)
  if (vine%npars()/=3) error stop 'C-vine truncation'
  print '(a)', 'test_cvine: PASS'
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
