! SPDX-License-Identifier: GPL-3.0-only
! Native Fortran translation of the computational core of rvinecopulib.
program example_cvine_fit
  use rvinecopulib
  implicit none
  type(cvine_model) :: truth,fit
  real(dp),allocatable :: sample(:,:)
  integer,parameter :: fams(4)=[bicop_indep,bicop_gaussian,bicop_clayton,bicop_frank]
  integer :: i,j
  allocate(sample(3,500))
  truth=make_cvine(3)
  truth%pair(1,2)=make_bicop(bicop_gaussian,0,[0.55_dp])
  truth%pair(1,3)=make_bicop(bicop_clayton,0,[1.0_dp])
  truth%pair(2,3)=make_bicop(bicop_frank,0,[2.0_dp])
  call seed_rng(777)
  call truth%simulate(size(sample,2),sample)
  fit=make_cvine(3)
  call fit%fit(sample,fams,criterion='bic')
  do i=1,2
    do j=i+1,3
      print '(a,i0,a,i0,a,a,a,i0,a,3f9.4)', 'pair(',i,',',j,'): ', &
        trim(fit%pair(i,j)%name()),' rot=',fit%pair(i,j)%rotation,' par=',fit%pair(i,j)%parameters
    end do
  end do
  print '(a,f12.4)', 'log likelihood: ',fit%loglik
end program
