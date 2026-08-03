! SPDX-License-Identifier: GPL-3.0-only
! Native Fortran translation of the computational core of rvinecopulib.
program test_dvine
  use rvinecopulib
  implicit none
  type(dvine_model) :: vine
  real(dp) :: x(3,1),z(3,1),back(3,1),sample(3,300),zr(3,300)
  integer :: i
  vine=make_dvine(3,[2,1,3])
  vine%pair(1,2)=make_bicop(bicop_gaussian,0,[0.5_dp])
  vine%pair(2,3)=make_bicop(bicop_clayton,0,[1.0_dp])
  vine%pair(1,3)=make_bicop(bicop_frank,0,[2.0_dp])
  x(:,1)=[0.2_dp,0.5_dp,0.8_dp]
  if (vine%pdf(x(:,1))<=0.0_dp) error stop 'D-vine density'
  call vine%rosenblatt(x,z)
  call vine%inverse_rosenblatt(z,back)
  if (maxval(abs(back-x))>2.0e-10_dp) error stop 'D-vine Rosenblatt roundtrip'
  call seed_rng(99)
  call vine%simulate(300,sample)
  call vine%rosenblatt(sample,zr)
  do i=1,3
    if (abs(sum(zr(i,:))/300.0_dp-0.5_dp)>0.07_dp) error stop 'D-vine transformed mean'
  end do
  print '(a)', 'test_dvine: PASS'
end program
