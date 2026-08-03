! SPDX-License-Identifier: GPL-3.0-only
program independence_example
  use iso_fortran_env, only : int64
  use mixedindtests
  implicit none
  real(dp) :: x(30,3)
  type(copula_test_result) :: test
  type(dependence_result) :: pairwise
  integer :: i

  do i = 1, 30
    x(i,1) = sin(0.7_dp*real(i,dp))
    x(i,2) = cos(1.1_dp*real(i,dp))
    x(i,3) = real(mod(7*i,13),dp)
  end do

  pairwise = EstDep(x)
  test = TestIndCopula(x, trunc_level=3, b=100, seed=2026_int64)

  print '(a,f10.4)', 'Pairwise LB p-value (%): ', pairwise%p_lb_rho
  print '(a,f10.4)', 'Combined CVM p-value (%): ', test%p_tn
  print '(a,f10.4)', 'Global Sn p-value (%): ', test%p_sn
end program independence_example
