! SPDX-License-Identifier: GPL-3.0-only
program demo_mixedindtests
  use iso_fortran_env, only : int64
  use mixedindtests
  implicit none
  real(dp), allocatable :: x(:),u(:)
  type(copula_test_result) :: test
  type(prepared_data_result) :: pd

  x=SimAR1Poisson([4.0_dp,0.25_dp],80,98765_int64)
  pd=preparedata(x)
  test=TestIndSerCopula(x,4,3,80,24680_int64)
  u=SimCopulaSeries('gaussian',12,0.4_dp,seed=13579_int64)

  print '(a,i0)', 'Unique values: ',size(pd%values)
  print '(a,f10.4)', 'Serial Sn: ',test%sn
  print '(a,f10.4)', 'Serial Sn p-value (%): ',test%p_sn
  print '(a,12f7.3)', 'Gaussian-copula uniforms: ',u
end program demo_mixedindtests
