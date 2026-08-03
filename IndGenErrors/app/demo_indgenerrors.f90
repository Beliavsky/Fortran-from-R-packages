! SPDX-License-Identifier: GPL-3.0-only
program demo_indgenerrors
  use indgenerrors
  implicit none
  integer, parameter :: n = 48
  real(dp) :: x(n), y(n), z(n)
  type(cvm_three_result) :: cvm
  type(four_lag_test_result) :: cor
  type(dependence_three_result) :: dep
  integer :: i

  do i = 1, n
    x(i) = sin(0.19_dp*real(i,dp))+0.02_dp*real(mod(i,5),dp)
    z(i) = sin(0.41_dp*real(i,dp))*cos(0.13_dp*real(i,dp))
  end do
  do i = 1, n
    y(i) = cos(0.29_dp*real(i,dp))+0.35_dp*x(modulo(i-2,n)+1)
  end do

  cvm = cvm_3series(x,y,z,3,1)
  cor = crosscor_3series(x,y,z,3,1)
  dep = crossdep_3series(x,y,z,3,1)

  print '(a,f10.4,2x,a,f8.5)', 'Combined CVM W =', cvm%wstat, 'p =', cvm%p_wstat
  print '(a,f10.4,2x,a,f8.5)', 'Combined correlation H =', cor%aggregate, &
    'p =', cor%p_aggregate
  print '(a,f10.4,2x,a,f8.5)', 'Combined Spearman H =', dep%spearman%aggregate, &
    'p =', dep%spearman%p_aggregate
end program demo_indgenerrors
