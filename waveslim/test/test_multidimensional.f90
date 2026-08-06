! SPDX-License-Identifier: BSD-3-Clause
! Modern Fortran computational translation of waveslim.
program test_multidimensional
  use waveslim
  use waveslim_test_support
  implicit none
  real(dp) :: image(32,32), cube(8,8,8)
  real(dp), allocatable :: image_back(:,:), cube_back(:,:,:), parts3(:,:,:,:)
  type(wavelet_transform_2d) :: w2
  type(wavelet_transform_3d) :: w3
  type(complex_wavelet_transform_2d) :: cd2
  integer :: i, j, k, level

  do j = 1, 32
    do i = 1, 32
      image(i,j) = sin(0.08_dp*real(i,dp)) + cos(0.13_dp*real(j,dp)) + &
        0.1_dp*sin(0.02_dp*real(i*j,dp))
    end do
  end do
  w2 = dwt_2d(image, 'haar', 3)
  image_back = idwt_2d(w2)
  call assert_close_scalar(maxval(abs(image-image_back)), 0.0_dp, &
    2.0e-12_dp, '2D DWT inverse')

  w2 = modwt_2d(image, 'la8', 3)
  image_back = imodwt_2d(w2)
  call assert_close_scalar(maxval(abs(image-image_back)), 0.0_dp, &
    2.0e-9_dp, '2D MODWT inverse')

  cd2 = dualtree_2d(image, 2)
  image_back = idualtree_2d(cd2)
  call assert_close_scalar(maxval(abs(image-image_back)), 0.0_dp, &
    8.0e-8_dp, '2D dual-tree inverse')

  do k = 1, 8
    do j = 1, 8
      do i = 1, 8
        cube(i,j,k) = sin(0.2_dp*real(i,dp)) + &
          cos(0.17_dp*real(j,dp)) + 0.05_dp*real(k,dp)
      end do
    end do
  end do
  w3 = dwt_3d(cube, 'haar', 2)
  cube_back = idwt_3d(w3)
  call assert_close_scalar(maxval(abs(cube-cube_back)), 0.0_dp, &
    5.0e-12_dp, '3D DWT inverse')

  w3 = modwt_3d(cube, 'haar', 2)
  cube_back = imodwt_3d(w3)
  call assert_close_scalar(maxval(abs(cube-cube_back)), 0.0_dp, &
    5.0e-12_dp, '3D MODWT inverse')

  parts3 = mra_3d(cube, 'haar', 2, 'modwt')
  cube_back = 0.0_dp
  do level = 1, size(parts3,4)
    cube_back = cube_back + parts3(:,:,:,level)
  end do
  call assert_close_scalar(maxval(abs(cube-cube_back)), 0.0_dp, &
    1.0e-11_dp, '3D MRA additive identity')

  write(*,'(a)') 'test_multidimensional: PASS'
end program test_multidimensional
