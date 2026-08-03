! SPDX-License-Identifier: GPL-3.0-only
! Native Fortran translation of the computational core of rvinecopulib.
program test_math_tools
  use rvinecopulib
  implicit none
  real(dp) :: p, x(2,5), u(2,5)
  integer :: i
  do i=1,9
    p=real(i,dp)/10.0_dp
    call assert_close(normal_cdf(normal_quantile(p)),p,2.0e-10_dp,'normal roundtrip')
    call assert_close(student_cdf(student_quantile(p,5.0_dp),5.0_dp),p,2.0e-10_dp,'student roundtrip')
  end do
  x(1,:)=[3.0_dp,1.0_dp,2.0_dp,2.0_dp,5.0_dp]
  x(2,:)=[5.0_dp,4.0_dp,3.0_dp,2.0_dp,1.0_dp]
  call pseudo_observations(x,u)
  call assert_close(u(1,1),4.0_dp/6.0_dp,1.0e-14_dp,'pseudo observation')
  call assert_close(u(1,3),2.5_dp/6.0_dp,1.0e-14_dp,'average tie rank')
  if (.not. is_valid_permutation([3,1,2])) error stop 'valid permutation failed'
  if (is_valid_permutation([1,1,3])) error stop 'invalid permutation accepted'
  print '(a)', 'test_math_tools: PASS'
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
