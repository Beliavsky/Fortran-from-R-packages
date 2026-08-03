! SPDX-License-Identifier: GPL-3.0-only
! Native Fortran translation of the computational core of rvinecopulib.
program test_vine_fit
  use rvinecopulib
  implicit none
  type(dvine_model) :: truth,fit
  real(dp),allocatable :: sample(:,:)
  integer,parameter :: fams(2)=[bicop_indep,bicop_gaussian]
  allocate(sample(3,450))
  truth=make_dvine(3)
  truth%pair(1,2)=make_bicop(bicop_gaussian,0,[0.6_dp])
  truth%pair(2,3)=make_bicop(bicop_gaussian,0,[-0.45_dp])
  truth%pair(1,3)=make_bicop(bicop_indep)
  call seed_rng(177)
  call truth%simulate(size(sample,2),sample)
  fit=make_dvine(3)
  call fit%fit(sample,fams,criterion='bic',allow_rotations=.false.)
  if (fit%pair(1,2)%family/=bicop_gaussian) error stop 'D-vine first edge selection'
  if (fit%pair(2,3)%family/=bicop_gaussian) error stop 'D-vine second edge selection'
  if (fit%loglik<=0.0_dp) error stop 'D-vine fit likelihood'
  if (fit%nobs/=size(sample,2)) error stop 'D-vine nobs'
  print '(a)', 'test_vine_fit: PASS'
end program
