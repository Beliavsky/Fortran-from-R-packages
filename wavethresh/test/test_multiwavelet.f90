! SPDX-License-Identifier: GPL-2.0-or-later
program test_multiwavelet
   use wavethresh
   implicit none
   real(dp), parameter :: tol = 2.0e-10_dp
   real(dp) :: x4(4)
   real(dp) :: x8(8)
   real(dp) :: x12(12)
   real(dp), allocatable :: coefficients(:,:)
   real(dp), allocatable :: reconstructed(:)
   real(dp), allocatable :: selected(:,:)
   type(mwd_t) :: object
   type(mwd_t) :: thresholded
   type(mwd_t) :: joint_object
   type(mwd_t) :: transformed
   logical :: ok
   integer :: i
   character(len=8), parameter :: geronimo_types(7) = [character(len=8) :: &
      "Minimal", "Identity", "Interp", "Xia", "Roach1", "Roach3", "default"]
   character(len=8), parameter :: donovan_types(4) = [character(len=8) :: &
      "Identity", "Linear", "Interp", "default"]

   do i = 1, 8
      x8(i) = real(i, dp) / 7.0_dp
   end do
   do i = 1, size(geronimo_types)
      coefficients = mprefilter(x8, trim(geronimo_types(i)), "Geronimo", 2, 20, 2, 2, 2)
      call assert_true(all(shape(coefficients) == [2, 20]), "Geronimo prefilter shape")
      reconstructed = mpostfilter(coefficients, trim(geronimo_types(i)), "Geronimo", 2, 2, 2, 2)
      call assert_vector_close(reconstructed, x8, tol, "Geronimo pre/postfilter round trip")
   end do

   transformed = mwd(x8, "Interp", "Geronimo")
   call assert_true(transformed%ok .and. transformed%nlevels == 2, "mwd Geronimo status")
   reconstructed = mwr(transformed)
   call assert_vector_close(reconstructed, x8, 5.0e-11_dp, "mwd/mwr Geronimo round trip")

   do i = 1, 4
      x4(i) = real(i, dp) / 5.0_dp
   end do
   transformed = mwd(x4, "Repeat", "Geronimo")
   call assert_true(transformed%ok .and. transformed%nlevels == 2, "mwd Geronimo Repeat status")
   reconstructed = mwr(transformed)
   call assert_vector_close(reconstructed, x4, 5.0e-11_dp, "mwd/mwr Geronimo Repeat round trip")

   coefficients = mprefilter(x4, "Repeat", "Geronimo", 2, 20, 2, 2, 2)
   reconstructed = mpostfilter(coefficients, "Repeat", "Geronimo", 2, 2, 2, 2)
   call assert_vector_close(reconstructed, x4, tol, "Geronimo Repeat round trip")

   do i = 1, 12
      x12(i) = sin(real(i, dp) / 5.0_dp)
   end do
   transformed = mwd(x12, "Interp", "Donovan3")
   call assert_true(transformed%ok .and. transformed%nlevels == 2, "mwd Donovan3 status")
   reconstructed = mwr(transformed)
   call assert_vector_close(reconstructed, x12, 5.0e-11_dp, "mwd/mwr Donovan3 round trip")

   do i = 1, size(donovan_types)
      coefficients = mprefilter(x12, trim(donovan_types(i)), "Donovan3", 2, 20, 3, 3, 2)
      call assert_true(all(shape(coefficients) == [3, 20]), "Donovan3 prefilter shape")
      reconstructed = mpostfilter(coefficients, trim(donovan_types(i)), "Donovan3", 3, 3, 2, 2)
      if (trim(donovan_types(i)) == "Linear") then
         call assert_vector_close(reconstructed, x12, 1.0e-7_dp, "Donovan3 Linear round trip")
      else
         call assert_vector_close(reconstructed, x12, tol, "Donovan3 pre/postfilter round trip")
      end if
   end do

   object%n_original = 8
   object%nlevels = 2
   object%filter = mfilter_select("Geronimo")
   allocate(object%detail(0:1))
   allocate(object%scaling(0:2))
   object%detail(0)%values = reshape([1.0_dp, 2.0_dp], [2, 1])
   object%detail(1)%values = reshape([3.0_dp, 4.0_dp, 5.0_dp, 6.0_dp], [2, 2])
   object%scaling(0)%values = reshape([7.0_dp, 8.0_dp], [2, 1])
   object%scaling(1)%values = reshape([9.0_dp, 10.0_dp, 11.0_dp, 12.0_dp], [2, 2])
   object%scaling(2)%values = reshape([(real(i, dp), i=1, 8)], [2, 4])
   object%ok = .true.
   object%message = "ok"

   selected = access_d_mwd(object, 1)
   call assert_true(all(shape(selected) == [2, 2]), "accessD.mwd shape")
   call assert_close(selected(2, 2), 6.0_dp, 0.0_dp, "accessD.mwd value")
   selected = access_c_mwd(object, 2)
   call assert_true(all(shape(selected) == [2, 4]), "accessC.mwd shape")
   call assert_close(selected(2, 4), 8.0_dp, 0.0_dp, "accessC.mwd value")

   call put_d_mwd(object, 1, reshape([13.0_dp, 14.0_dp, 15.0_dp, 16.0_dp], [2, 2]), ok)
   call assert_true(ok, "putD.mwd status")
   selected = access_d_mwd(object, 1)
   call assert_close(selected(1, 2), 15.0_dp, 0.0_dp, "putD.mwd value")

   thresholded = threshold_mwd(object, levels=[1], threshold_type="hard", policy="manual", &
      value=[14.5_dp, 15.5_dp], bivariate=.false.)
   call assert_true(thresholded%ok, "threshold.mwd componentwise status")
   selected = access_d_mwd(thresholded, 1)
   call assert_close(selected(1, 1), 0.0_dp, 0.0_dp, "threshold.mwd component 1 zero")
   call assert_close(selected(1, 2), 15.0_dp, 0.0_dp, "threshold.mwd component 1 keep")
   call assert_close(selected(2, 1), 0.0_dp, 0.0_dp, "threshold.mwd component 2 zero")
   call assert_close(selected(2, 2), 16.0_dp, 0.0_dp, "threshold.mwd component 2 keep")

   joint_object = object
   joint_object%detail(1)%values = reshape([1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, &
      -1.0_dp, 0.0_dp, 0.0_dp, -1.0_dp], [2, 4])
   thresholded = threshold_mwd(joint_object, levels=[1], threshold_type="hard", policy="manual", &
      value=[2.0_dp], robust=.false., bivariate=.true.)
   selected = access_d_mwd(thresholded, 1)
   call assert_close(maxval(abs(selected)), 0.0_dp, 0.0_dp, "threshold.mwd joint chi-square")

   call put_c_mwd(object, 0, reshape([17.0_dp, 18.0_dp], [2, 1]), ok)
   call assert_true(ok, "putC.mwd status")
   selected = access_c_mwd(object, 0)
   call assert_close(selected(2, 1), 18.0_dp, 0.0_dp, "putC.mwd value")
   call assert_true(nlevels_mwd(object) == 2, "nlevelsWT.mwd typed depth")

   print *, "test_multiwavelet: PASS"

contains

   subroutine assert_true(condition, message)
      logical, intent(in) :: condition !! Condition that must evaluate true for the test to pass.
      character(len=*), intent(in) :: message !! Diagnostic label printed when the assertion fails.
      if (.not. condition) then
         print *, "FAIL: ", trim(message)
         error stop 1
      end if
   end subroutine assert_true

   subroutine assert_close(actual, expected, tolerance, message)
      real(dp), intent(in) :: actual !! Scalar value produced by the translated implementation.
      real(dp), intent(in) :: expected !! Deterministic scalar reference value.
      real(dp), intent(in) :: tolerance !! Maximum permitted absolute difference from the reference.
      character(len=*), intent(in) :: message !! Diagnostic label printed when the assertion fails.
      if (abs(actual - expected) > tolerance) then
         print *, "FAIL: ", trim(message), actual, expected
         error stop 1
      end if
   end subroutine assert_close

   subroutine assert_vector_close(actual, expected, tolerance, message)
      real(dp), intent(in) :: actual(:) !! Vector produced by the translated implementation.
      real(dp), intent(in) :: expected(:) !! Deterministic vector of reference values.
      real(dp), intent(in) :: tolerance !! Maximum permitted elementwise absolute error.
      character(len=*), intent(in) :: message !! Diagnostic label printed when the assertion fails.
      if (size(actual) /= size(expected)) then
         print *, "FAIL size: ", trim(message)
         error stop 1
      end if
      if (size(actual) > 0) then
         if (maxval(abs(actual - expected)) > tolerance) then
            print *, "FAIL: ", trim(message)
            error stop 1
         end if
      end if
   end subroutine assert_vector_close

end program test_multiwavelet
