! SPDX-License-Identifier: GPL-2.0-or-later
program test_uniroot_all
  use rootsolve, only : dp, uniroot_all
  implicit none
  real(dp),allocatable::r(:)
  real(dp),parameter::pi=acos(-1.0_dp)
  integer::i
  r=uniroot_all(fun,0.0_dp,10.0_dp,ngrid=200,tol=1e-10_dp)
  if(size(r)/=6) error stop 1
  do i=1,6
    if(abs(r(i)-(pi/4.0_dp+real(i-1,dp)*pi/2.0_dp))>1e-8_dp) error stop 2
  end do
  print *, 'test_uniroot_all: PASS'
contains
  real(dp) function fun(x) result(v)
    real(dp),intent(in)::x
    v=cos(2.0_dp*x)**3
  end function fun
end program test_uniroot_all
