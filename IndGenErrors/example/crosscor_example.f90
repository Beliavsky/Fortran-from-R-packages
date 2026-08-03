! SPDX-License-Identifier: GPL-3.0-only
program crosscor_example
  use indgenerrors
  implicit none
  integer, parameter :: n = 60
  real(dp) :: x(n), y(n)
  type(lag_test_result) :: result
  integer :: i

  do i = 1, n
    x(i) = sin(0.17_dp*real(i,dp))
  end do
  do i = 1, n
    y(i) = 0.65_dp*x(modulo(i-3,n)+1)+cos(0.43_dp*real(i,dp))
  end do
  result = crosscor_2series(x,y,5)
  print '(a,*(1x,i3))', 'lags:', result%lags(1,:)
  print '(a,*(1x,f8.4))', 'correlations:', result%stat
  print '(a,f10.4,2x,a,f8.5)', 'Wald statistic =', result%aggregate, &
    'p =', result%p_aggregate
end program crosscor_example
