! SPDX-License-Identifier: GPL-3.0-only
program cvm_example
  use indgenerrors
  implicit none
  integer, parameter :: n = 40
  real(dp) :: x(n), y(n)
  type(cvm_test_result) :: result
  integer :: i

  do i = 1, n
    x(i) = sin(0.31_dp*real(i,dp))+0.05_dp*real(mod(i,5),dp)
    y(i) = cos(0.47_dp*real(i,dp))+0.04_dp*real(mod(3*i,7),dp)
  end do
  result = cvm_2series(x,y,3)
  print '(a,*(1x,f9.5))', 'CVM statistics:', result%cvm
  print '(a,f10.5,2x,a,f8.5)', 'W =', result%wstat, 'p =', result%p_wstat
  print '(a,f10.5,2x,a,f8.5)', 'F =', result%fstat, 'p =', result%p_fstat
end program cvm_example
