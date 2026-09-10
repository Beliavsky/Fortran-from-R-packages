! SPDX-License-Identifier: GPL-2.0-or-later
program test_nd
   use wavethresh
   implicit none
   real(dp), parameter :: tol = 2.0e-11_dp
   real(dp) :: image(4, 4)
   real(dp) :: cube(4, 4, 4)
   real(dp), allocatable :: image_back(:,:)
   real(dp), allocatable :: cube_back(:,:,:)
   type(imwd_t) :: image_wd
   type(imwd_t) :: image_wst
   type(imwd_t) :: compressed
   type(imwd_t) :: restored
   type(imwd_t) :: thresholded
   type(wd3d_t) :: cube_wd
   integer :: i
   integer :: j
   integer :: k

   do j = 1, 4
      do i = 1, 4
         image(i, j) = real(i + 4 * (j - 1), dp) / 7.0_dp
      end do
   end do
   image_wd = imwd(image, filter_number=1.0_dp, family="DaubExPhase", nlevels=2)
   call assert_true(image_wd%ok, "imwd status")
   image_back = imwr(image_wd)
   call assert_close(maxval(abs(image_back - image)), 0.0_dp, tol, "imwd round trip")
   compressed = compress_imwd(image_wd)
   restored = uncompress_imwd(compressed)
   image_back = imwr(restored)
   call assert_close(maxval(abs(image_back - image)), 0.0_dp, tol, "imwd compressed storage round trip")
   thresholded = threshold_imwd(compressed, policy="manual", value=0.0_dp)
   image_back = imwr(thresholded)
   call assert_close(maxval(abs(image_back - image)), 0.0_dp, tol, "threshold.imwdc zero threshold")

   image_wst = wst2d(image, filter_number=1.0_dp, family="DaubExPhase", nlevels=2)
   call assert_true(image_wst%ok, "wst2D status")
   image_back = iwst2d(image_wst)
   call assert_close(maxval(abs(image_back - image)), 0.0_dp, 5.0e-11_dp, "wst2D round trip")

   do k = 1, 4
      do j = 1, 4
         do i = 1, 4
            cube(i, j, k) = sin(real(i + 2 * j + 3 * k, dp))
         end do
      end do
   end do
   cube_wd = wd3d(cube, filter_number=1.0_dp, family="DaubExPhase", nlevels=2)
   call assert_true(cube_wd%ok, "wd3D status")
   cube_back = wr3d(cube_wd)
   call assert_close(maxval(abs(cube_back - cube)), 0.0_dp, 8.0e-11_dp, "wd3D round trip")

   print *, "test_nd: PASS"
contains
   subroutine assert_true(condition, message)
      logical, intent(in) :: condition !! Condition that must hold.
      character(len=*), intent(in) :: message !! Diagnostic label for a failed assertion.
      if (.not. condition) then
         print *, "FAIL: ", trim(message)
         error stop 1
      end if
   end subroutine assert_true

   subroutine assert_close(actual, expected, tolerance, message)
      real(dp), intent(in) :: actual !! Value produced by the translation.
      real(dp), intent(in) :: expected !! Deterministic reference value.
      real(dp), intent(in) :: tolerance !! Maximum accepted absolute error.
      character(len=*), intent(in) :: message !! Diagnostic label for a failed assertion.
      if (abs(actual - expected) > tolerance) then
         print *, "FAIL: ", trim(message), actual, expected
         error stop 1
      end if
   end subroutine assert_close
end program test_nd
