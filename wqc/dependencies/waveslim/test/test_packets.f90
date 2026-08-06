! SPDX-License-Identifier: BSD-3-Clause
! Modern Fortran computational translation of waveslim.
program test_packets
  use waveslim
  use waveslim_test_support
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  implicit none
  real(dp) :: x(128)
  real(dp), allocatable :: y(:), entropy(:)
  logical, allocatable :: css(:), cpgram(:), port(:), basis(:)
  type(packet_transform) :: tree, mtree
  integer :: i

  do i = 1, size(x)
    x(i) = sin(0.09_dp*real(i,dp)) + 0.2_dp*cos(0.31_dp*real(i,dp))
  end do
  tree = dwpt(x, 'la8', 4)
  call assert_true(tree%status%ok(), 'DWPT status')
  y = idwpt(tree)
  call assert_close_array(y, x, 1.0e-9_dp, 'DWPT inverse')

  mtree = modwpt(x, 'la8', 3)
  call assert_true(mtree%status%ok(), 'MODWPT status')
  call assert_true(size(mtree%level(3)%node) == 8, 'MODWPT node count')
  call assert_true(size(mtree%level(3)%node(1)%values) == size(x), &
    'MODWPT undecimated node length')

  entropy = entropy_test(tree)
  css = css_test(tree)
  cpgram = cpgram_test(tree)
  port = portmanteau_test(tree)
  call assert_true(size(entropy) == 30, 'packet entropy result size')
  call assert_true(size(css) == 30 .and. size(cpgram) == 30, &
    'packet white-noise test sizes')
  call assert_true(size(port) == 30, 'packet portmanteau size')
  call assert_true(all(ieee_is_finite(entropy)), 'packet entropy finite')

  basis = packet_basis(tree, [1,2,2], [0,2,3])
  call assert_true(count(basis) == 3, 'packet basis selection')

  write(*,'(a)') 'test_packets: PASS'
end program test_packets
