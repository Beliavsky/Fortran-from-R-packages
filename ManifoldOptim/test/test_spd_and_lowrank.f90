! SPDX-License-Identifier: GPL-2.0-or-later
! Modern Fortran translation of the ManifoldOptim/ROPTLIB computational interface.
! See LICENSE and UPSTREAM_PROVENANCE.md for upstream copyright and provenance.
program test_spd_and_lowrank
  use manifoldoptim
  implicit none
  type(manifold_domain)::spd,lr
  real(dp),allocatable::x(:),eta(:),eta_proj(:),y(:)
  integer::n,m,r,i
  logical::ok
  n=3
  allocate(spd%component(1))
  spd%component(1)=make_component(MANI_SPD,n)
  allocate(x(n*n),eta(n*n),eta_proj(n*n),y(n*n))
  x=0.0_dp
  do i=1,n
  x(i+(i-1)*n)=1.0_dp
  end do
  eta=[0.2_dp,0.1_dp,0.0_dp, 0.1_dp,-0.1_dp,0.05_dp, 0.0_dp,0.05_dp,0.3_dp]
  call project_tangent(spd,x,eta,eta_proj)
  call retract_point(spd,x,eta_proj,y,ok)
  if(.not.ok.or..not.point_is_valid(spd,y,1.0e-10_dp))error stop 'SPD retraction'
  deallocate(x,eta,eta_proj,y)
  n=5
  m=4
  r=2
  allocate(lr%component(1))
  lr%component(1)=make_component(MANI_LOWRANK,n,m=m,p=r)
  allocate(x(lr%length()),eta(lr%length()),eta_proj(lr%length()),y(lr%length()))
  call random_manifold_point(lr,x)
  call random_number(eta)
  eta=0.05_dp*(2.0_dp*eta-1.0_dp)
  call project_tangent(lr,x,eta,eta_proj)
  call retract_point(lr,x,eta_proj,y,ok)
  if(.not.ok.or..not.point_is_valid(lr,y,1.0e-8_dp))error stop 'lowrank retraction'
  write(*,*)'PASS test_spd_and_lowrank'
end program
