! SPDX-License-Identifier: GPL-3.0-only
! Native Fortran translation of the computational core of rvinecopulib.
program demo_rvinecopulib
  use rvinecopulib
  implicit none
  type(dvine_model) :: truth,fit
  type(bicop_model) :: pair_fit
  real(dp),allocatable :: raw(:,:),u(:,:),sample(:,:),z(:,:)
  integer,parameter :: fams(5)=[bicop_indep,bicop_gaussian,bicop_student, &
                                bicop_clayton,bicop_frank]
  integer :: i
  allocate(raw(3,400),u(3,400),sample(3,400),z(3,400))
  call seed_rng(8128)
  truth=make_dvine(3)
  truth%pair(1,2)=make_bicop(bicop_gaussian,0,[0.65_dp])
  truth%pair(2,3)=make_bicop(bicop_clayton,0,[1.1_dp])
  truth%pair(1,3)=make_bicop(bicop_frank,0,[2.5_dp])
  call truth%simulate(400,u)
  do i=1,400
    raw(1,i)=normal_quantile(u(1,i))
    raw(2,i)=2.0_dp+0.5_dp*normal_quantile(u(2,i))
    raw(3,i)=exp(normal_quantile(u(3,i)))
  end do
  call pseudo_observations(raw,u)
  call select_bicop(u(1:2,:),pair_fit,fams,criterion='bic')
  print '(a,a,a,i0)', 'bivariate selection: ',trim(pair_fit%name()), &
    ', rotation ',pair_fit%rotation
  fit=make_dvine(3)
  call fit%fit(u,fams,criterion='bic')
  print '(a,f12.4,a,f12.4)', 'D-vine logLik = ',fit%loglik,', BIC = ',fit%bic()
  call fit%simulate(400,sample)
  call fit%rosenblatt(sample,z)
  print '(a,3f9.5)', 'mean Rosenblatt coordinates: ', &
    sum(z(1,:))/400.0_dp,sum(z(2,:))/400.0_dp,sum(z(3,:))/400.0_dp
end program
