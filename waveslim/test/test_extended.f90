! SPDX-License-Identifier: BSD-3-Clause
! Modern Fortran computational translation of waveslim.
program test_extended
  use waveslim
  use waveslim_test_support
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  implicit none
  real(dp) :: image(16,16), x(256), variance_signal(512)
  real(dp), allocatable :: reconstructed(:,:), boot1(:), boot2(:)
  type(packet_transform_2d) :: tree2d
  type(variance_change_result) :: changes
  integer :: i, j

  do j = 1, size(image,2)
    do i = 1, size(image,1)
      image(i,j) = sin(0.17_dp*real(i,dp)) + cos(0.11_dp*real(j,dp)) &
        + 0.01_dp*real(i*j,dp)
    end do
  end do
  tree2d = dwpt_2d(image, 'haar', 2)
  call assert_true(tree2d%status%ok(), '2D DWPT status')
  call assert_true(tree2d%levels() == 2, '2D DWPT level count')
  call assert_true(size(tree2d%level(2)%node) == 16, '2D DWPT node count')
  reconstructed = idwpt_2d(tree2d)
  call assert_true(all(shape(reconstructed) == shape(image)), &
    '2D DWPT inverse shape')
  call assert_true(maxval(abs(reconstructed-image)) < 1.0e-10_dp, &
    '2D DWPT inverse')

  do i = 1, size(x)
    x(i) = sin(0.13_dp*real(i,dp)) + 0.4_dp*cos(0.37_dp*real(i,dp))
  end do
  boot1 = dwpt_boot(x, 'haar', 4, 4321_i8)
  boot2 = dwpt_boot(x, 'haar', 4, 4321_i8)
  call assert_true(size(boot1) == size(x), 'DWPT bootstrap length')
  call assert_true(all(ieee_is_finite(boot1)), 'DWPT bootstrap finite')
  call assert_close_array(boot1, boot2, 0.0_dp, 'DWPT bootstrap seeded')

  do i = 1, size(variance_signal)
    if (i <= size(variance_signal)/2) then
      variance_signal(i) = sin(0.29_dp*real(i,dp)) &
        + 0.3_dp*cos(0.47_dp*real(i,dp))
    else
      variance_signal(i) = 4.0_dp*(sin(0.29_dp*real(i,dp)) &
        + 0.3_dp*cos(0.47_dp*real(i,dp)))
    end if
  end do
  changes = testing_hov(variance_signal, 'haar', 3, 16)
  call assert_true(changes%status%ok(), 'variance-change status')
  call assert_true(changes%count() > 0, 'variance change detected')
  call assert_true(all(changes%level >= 1 .and. changes%level <= 3), &
    'variance-change levels')
  call assert_true(all(changes%statistic > 0.0_dp), &
    'variance-change statistics')

  write(*,'(a)') 'test_extended: PASS'
end program test_extended
