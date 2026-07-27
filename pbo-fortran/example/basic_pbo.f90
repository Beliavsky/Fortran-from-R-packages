! SPDX-License-Identifier: MIT
! Copyright (c) 2014 Matthew R. Barry
program basic_pbo
  use pbo, only : dp, pbo_result, selection_result, compute_pbo, column_mean, &
    selection_frequencies
  implicit none
  type(pbo_result) :: result
  type(selection_result) :: frequencies
  real(dp) :: x(8,3)
  integer :: i

  do i = 1, 8
    x(i,1) = 0.1_dp * real(i,dp)
    x(i,2) = 0.9_dp - 0.08_dp * real(i,dp)
    x(i,3) = 0.2_dp * sin(real(i,dp))
  end do
  call compute_pbo(x, 4, column_mean, result, threshold=0.0_dp)
  if (.not. result%success) error stop result%message
  call selection_frequencies(result, frequencies)
  print '(a,f8.4)', 'PBO = ', result%phi
  print '(a,*(i0,1x))', 'Selection order: ', frequencies%strategy
end program basic_pbo
