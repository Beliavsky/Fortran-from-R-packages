! SPDX-License-Identifier: GPL-2.0-or-later
program test_threshold_extra
   use wavethresh
   implicit none
   real(dp), parameter :: tol = 5.0e-12_dp
   real(dp) :: data(8)
   real(dp) :: outlier(8)
   real(dp), allocatable :: cc(:)
   real(dp), allocatable :: thresholds(:)
   real(dp) :: threshold
   type(wd_t) :: object
   type(wd_t) :: thresholded
   integer :: level

   data = [0.1_dp, -0.3_dp, 0.8_dp, -1.2_dp, 2.0_dp, -0.7_dp, 1.5_dp, 0.05_dp]
   threshold = to_getthrda1(data, 0.05_dp)
   call assert_close(threshold, 2.0_dp, tol, "TOgetthrda1")
   cc = to_one_by_one1(sort_squares(data), 0.05_dp)
   call assert_true(size(cc) == 1, "TOonebyone1 sequence length")

   outlier = [0.1_dp, 0.1_dp, 0.1_dp, 0.1_dp, 0.1_dp, 0.1_dp, 0.1_dp, 5.0_dp]
   cc = to_one_by_one2(outlier**2, 0.05_dp)
   call assert_true(size(cc) == 2, "TOonebyone2 sequence length")
   threshold = to_getthrda2(outlier**2, 0.05_dp)
   call assert_close(threshold, 0.1_dp, tol, "TOgetthrda2")

   object = wd(data, filter_number=1.0_dp, family="DaubExPhase")
   call assert_true(object%ok, "TOthresh input transform")
   thresholds = to_threshda1_thresholds(object, 0.05_dp)
   call assert_true(size(thresholds) == object%nlevels, "TOthreshda1 threshold count")
   thresholded = to_threshda1(object, 0.05_dp)
   call assert_true(thresholded%ok, "TOthreshda1 status")
   do level = 0, object%nlevels - 1
      call assert_true(all(abs(thresholded%detail(level)%values) <= abs(object%detail(level)%values) + tol), &
         "TOthreshda1 shrinkage")
   end do

   thresholds = to_threshda2_thresholds(object, 0.05_dp)
   call assert_true(size(thresholds) == object%nlevels, "TOthreshda2 threshold count")
   thresholded = to_threshda2(object, 0.05_dp)
   call assert_true(thresholded%ok, "TOthreshda2 status")
   do level = 0, object%nlevels - 1
      call assert_true(all(abs(thresholded%detail(level)%values) <= abs(object%detail(level)%values) + tol), &
         "TOthreshda2 shrinkage")
   end do

   print *, "test_threshold_extra: PASS"

contains

   function sort_squares(values) result(sorted)
      real(dp), intent(in) :: values(:) !! Coefficients whose squared values are returned in ascending order.
      real(dp), allocatable :: sorted(:)
      real(dp) :: key
      integer :: i
      integer :: j

      sorted = values**2
      do i = 2, size(sorted)
         key = sorted(i)
         j = i - 1
         do while (j >= 1)
            if (sorted(j) <= key) exit
            sorted(j + 1) = sorted(j)
            j = j - 1
         end do
         sorted(j + 1) = key
      end do
   end function sort_squares

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

end program test_threshold_extra
