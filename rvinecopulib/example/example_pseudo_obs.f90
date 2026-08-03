! SPDX-License-Identifier: GPL-3.0-only
! Native Fortran translation of the computational core of rvinecopulib.
program example_pseudo_obs
  use rvinecopulib
  implicit none
  real(dp) :: x(2,6),u(2,6)
  integer :: i
  x(1,:)=[10.0_dp,8.0_dp,11.0_dp,9.0_dp,9.0_dp,14.0_dp]
  x(2,:)=[2.0_dp,5.0_dp,4.0_dp,1.0_dp,3.0_dp,6.0_dp]
  call pseudo_observations(x,u)
  print '(a)', 'column:       1         2         3         4         5         6'
  do i=1,2
    print '(a,i0,a,6f10.5)', 'margin ',i,': ',u(i,:)
  end do
end program
