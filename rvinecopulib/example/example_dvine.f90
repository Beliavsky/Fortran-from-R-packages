! SPDX-License-Identifier: GPL-3.0-only
! Native Fortran translation of the computational core of rvinecopulib.
program example_dvine
  use rvinecopulib
  implicit none
  type(dvine_model) :: vine
  real(dp) :: sample(3,5),z(3,5)
  integer :: j
  vine=make_dvine(3)
  vine%pair(1,2)=make_bicop(bicop_gaussian,0,[0.6_dp])
  vine%pair(2,3)=make_bicop(bicop_clayton,0,[1.2_dp])
  vine%pair(1,3)=make_bicop(bicop_frank,0,[2.0_dp])
  call seed_rng(2026)
  call vine%simulate(5,sample)
  call vine%rosenblatt(sample,z)
  do j=1,5
    print '(a,i0,a,3f9.5,a,3f9.5)', 'obs ',j,': ',sample(:,j),'  ->  ',z(:,j)
  end do
end program
