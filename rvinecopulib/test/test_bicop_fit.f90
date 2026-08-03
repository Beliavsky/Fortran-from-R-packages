! SPDX-License-Identifier: GPL-3.0-only
! Native Fortran translation of the computational core of rvinecopulib.
program test_bicop_fit
  use rvinecopulib
  implicit none
  type(bicop_model) :: truth,fit,selected
  real(dp),allocatable :: sample(:,:)
  allocate(sample(2,500))
  call seed_rng(731)
  truth=make_bicop(bicop_gaussian,0,[0.65_dp])
  call truth%simulate(size(sample,2),sample)
  call fit_bicop(sample,bicop_gaussian,0,fit,max_iter=400)
  if (abs(fit%parameters(1)-0.65_dp)>0.10_dp) then
    print *,fit%parameters(1)
    error stop 'Gaussian MLE inaccurate'
  end if
  call select_bicop(sample,selected,[bicop_indep,bicop_gaussian,bicop_clayton],criterion='bic')
  if (selected%family/=bicop_gaussian) then
    print *,selected%family,selected%rotation,selected%parameters
    error stop 'family selection failed'
  end if
  if (selected%loglik<=0.0_dp) error stop 'loglik not retained'
  print '(a)', 'test_bicop_fit: PASS'
end program
