! SPDX-License-Identifier: GPL-3.0-only
! Native Fortran translation of the computational core of rvinecopulib.
program example_bicop_fit
  use rvinecopulib
  implicit none
  type(bicop_model) :: truth,fit
  real(dp),allocatable :: u(:,:)
  allocate(u(2,600))
  call seed_rng(1234)
  truth=make_bicop(bicop_gaussian,0,[0.70_dp])
  call truth%simulate(size(u,2),u)
  call select_bicop(u,fit,[bicop_indep,bicop_gaussian,bicop_student, &
                           bicop_clayton,bicop_gumbel],criterion='bic')
  print '(a,a)', 'selected family: ',trim(fit%name())
  print '(a,i0)', 'rotation:        ',fit%rotation
  print '(a,3f10.5)', 'parameters:      ',fit%parameters
  print '(a,f12.4)', 'log likelihood:  ',fit%loglik
  print '(a,f12.4)', 'BIC:             ',fit%bic()
end program
