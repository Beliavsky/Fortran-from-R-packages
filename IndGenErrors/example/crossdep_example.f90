! SPDX-License-Identifier: GPL-3.0-only
program crossdep_example
  use indgenerrors
  implicit none
  integer, parameter :: n = 50
  real(dp) :: x(n), y(n)
  type(dependence_two_result) :: result
  integer :: i

  do i = 1, n
    x(i) = sin(0.21_dp*real(i,dp))+0.01_dp*real(i,dp)
    y(i) = sign(1.0_dp,x(i))*cos(0.37_dp*real(i,dp))**2+0.02_dp*real(mod(i,4),dp)
  end do
  result = crossdep_2series(x,y,3)
  print '(a,*(1x,f8.4))', 'Spearman:', result%spearman%stat
  print '(a,*(1x,f8.4))', 'van der Waerden:', result%vdw%stat
  print '(a,*(1x,f8.4))', 'Savage:', result%savage%stat
  print '(a,f10.4,2x,a,f8.5)', 'Spearman H =', result%spearman%aggregate, &
    'p =', result%spearman%p_aggregate
end program crossdep_example
