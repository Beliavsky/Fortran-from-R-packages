! SPDX-License-Identifier: BSD-3-Clause
! Modern Fortran computational translation of waveslim.
program image_packet_example
  use waveslim
  implicit none
  real(dp) :: image(16,16)
  real(dp), allocatable :: reconstructed(:,:)
  type(packet_transform_2d) :: tree
  integer :: i, j

  do j = 1, size(image,2)
    do i = 1, size(image,1)
      image(i,j) = sin(0.2_dp*real(i,dp)) + cos(0.15_dp*real(j,dp))
    end do
  end do

  tree = dwpt_2d(image, 'haar', 2)
  if (.not. tree%status%ok()) error stop trim(tree%status%message)
  reconstructed = idwpt_2d(tree)

  write(*,'(a,i0)') 'Packet levels: ', tree%levels()
  write(*,'(a,i0)') 'Nodes at level 2: ', size(tree%level(2)%node)
  write(*,'(a,es12.4)') 'Maximum reconstruction error: ', &
    maxval(abs(reconstructed-image))
end program image_packet_example
