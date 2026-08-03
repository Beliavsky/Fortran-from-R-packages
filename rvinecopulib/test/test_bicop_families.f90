! SPDX-License-Identifier: GPL-3.0-only
! Native Fortran translation of the computational core of rvinecopulib.
program test_bicop_families
  use rvinecopulib
  implicit none
  type(bicop_model) :: m
  integer :: f
  real(dp) :: par(3),c,d,h,tau
  do f=bicop_joe,bicop_tawn
    select case(f)
    case(bicop_joe); par=[1.5_dp,0.0_dp,0.0_dp]
    case(bicop_bb1); par=[1.0_dp,1.5_dp,0.0_dp]
    case(bicop_bb6); par=[1.5_dp,1.2_dp,0.0_dp]
    case(bicop_bb7); par=[1.5_dp,1.0_dp,0.0_dp]
    case(bicop_bb8); par=[1.5_dp,0.8_dp,0.0_dp]
    case(bicop_tawn); par=[0.8_dp,0.9_dp,1.5_dp]
    end select
    m=make_bicop(f,0,par(1:family_parameter_count(f)))
    c=m%cdf(0.3_dp,0.7_dp)
    d=m%pdf(0.3_dp,0.7_dp)
    h=m%hfunc1(0.3_dp,0.7_dp)
    if (c<0.0_dp .or. c>0.3_dp) error stop 'cdf outside Frechet bounds'
    if (.not.(d>0.0_dp)) error stop 'nonpositive density'
    if (h<=0.0_dp .or. h>=1.0_dp) error stop 'invalid h function'
    tau=m%tau()
    if (tau<0.0_dp .or. tau>1.0_dp) error stop 'invalid tau'
  end do
  m=make_bicop(bicop_student,0,[0.4_dp,6.0_dp])
  if (m%pdf(0.2_dp,0.8_dp)<=0.0_dp) error stop 'student density'
  call assert_close(m%hinv1(0.2_dp,m%hfunc1(0.2_dp,0.8_dp)),0.8_dp,3.0e-10_dp,'student inverse')
  print '(a)', 'test_bicop_families: PASS'
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
