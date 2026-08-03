! SPDX-License-Identifier: GPL-3.0-only
program serial_randomness_example
  use iso_fortran_env, only : int64
  use mixedindtests
  implicit none
  real(dp), allocatable :: x(:)
  type(serial_dependence_result) :: pairwise
  type(copula_test_result) :: test

  x = SimAR1Poisson([3.0_dp,0.35_dp],120,1234_int64)
  pairwise = EstDepSerial(x,5)
  test = TestIndSerCopula(x,4,trunc_level=3,b=100,seed=4321_int64)

  print '(a,5f9.4)', 'Serial rho: ', pairwise%rho
  print '(a,f10.4)', 'Serial LB p-value (%): ', pairwise%p_lb_rho
  print '(a,f10.4)', 'Copula combined p-value (%): ', test%p_tn
end program serial_randomness_example
