! SPDX-License-Identifier: GPL-3.0-only
! Native Fortran translation of the computational core of rvinecopulib.
program example_bicop_evaluate
  use rvinecopulib
  implicit none
  type(bicop_model) :: cop
  real(dp) :: u,v
  cop=make_bicop(bicop_clayton,0,[2.0_dp])
  u=0.35_dp; v=0.70_dp
  print '(a,a)', 'family: ',trim(cop%name())
  print '(a,f10.6)', 'cdf:    ',cop%cdf(u,v)
  print '(a,f10.6)', 'pdf:    ',cop%pdf(u,v)
  print '(a,f10.6)', 'h1:     ',cop%hfunc1(u,v)
  print '(a,f10.6)', 'h2:     ',cop%hfunc2(u,v)
  print '(a,f10.6)', 'tau:    ',cop%tau()
end program
