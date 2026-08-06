! SPDX-License-Identifier: BSD-3-Clause
! Modern Fortran computational translation of waveslim.
program test_errors
  use waveslim
  use waveslim_test_support
  implicit none
  real(dp) :: x(30)
  type(wavelet_transform) :: wt
  type(packet_transform) :: tree
  type(wavelet_filter_type) :: filter
  type(status_type) :: status

  x = 0.0_dp
  wt = dwt(x, 'haar', 4)
  call assert_true(.not. wt%status%ok(), 'invalid dyadic length rejected')
  wt = modwt(x, 'haar', 6)
  call assert_true(.not. wt%status%ok(), 'excessive MODWT depth rejected')
  tree = dwpt(x, 'haar', 4)
  call assert_true(.not. tree%status%ok(), 'invalid packet depth rejected')
  filter = wave_filter('not-a-filter', status)
  call assert_true(.not. status%ok(), 'unknown filter rejected')
  call assert_true(filter%length() == 0, 'unknown filter is empty')

  write(*,'(a)') 'test_errors: PASS'
end program test_errors
