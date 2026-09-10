program test_irregular
   use wavethresh, only : accessc, dp, grid_data_t, irregular_wavelet_t, irregwd, makegrid, threshold_irregwd
   implicit none

   type(grid_data_t) :: grid
   type(irregular_wavelet_t) :: object
   type(irregular_wavelet_t) :: thresholded
   real(dp) :: locations(8)
   real(dp) :: observations(8)
   real(dp), allocatable :: variances(:)
   real(dp), allocatable :: expected(:)
   real(dp) :: center
   real(dp) :: noise
   real(dp) :: cutoff
   integer :: i

   do i = 1, 8
      locations(i) = (real(i, dp) - 0.5_dp) / 8.0_dp
      observations(i) = sin(0.7_dp * real(i, dp)) + 0.15_dp * cos(0.23_dp * real(i * i, dp))
   end do
   grid = makegrid(locations, observations, 8)
   call check(grid%ok, "makegrid status")
   object = irregwd(grid, 1.0_dp, "DaubExPhase", "periodic")
   call check(object%ok, "irregwd status")
   do i = 0, object%transform%nlevels - 1
      variances = accessc(object, i)
      call check(size(variances) == size(object%transform%detail(i)%values), "accessc level shape")
      call check_close(maxval(abs(variances - 1.0_dp)), 0.0_dp, 3.0e-12_dp, "identity-grid variance")
   end do

   expected = object%transform%detail(2)%values
   center = sum(expected) / real(size(expected), dp)
   noise = sqrt(sum((expected - center)**2) / real(size(expected) - 1, dp))
   cutoff = 0.5_dp * noise
   where (abs(expected) <= cutoff) expected = 0.0_dp
   thresholded = threshold_irregwd(object, [2], "hard", "manual", 0.5_dp)
   call check(thresholded%ok, "threshold.irregwd status")
   call check_close(maxval(abs(thresholded%transform%detail(2)%values - expected)), 0.0_dp, 3.0e-12_dp, &
      "variance-adjusted manual threshold")

   thresholded = threshold_irregwd(object, [2], "soft", "manual", 0.0_dp)
   call check_close(maxval(abs(thresholded%transform%detail(2)%values - &
      object%transform%detail(2)%values)), 0.0_dp, 0.0_dp, "zero soft threshold")
   variances = accessc(object, -1)
   call check(size(variances) == 0, "accessc rejects invalid level")
   print *, "test_irregular: PASS"

contains

   subroutine check(condition, label)
      !! Stops the test program when a logical assertion fails.
      logical, intent(in) :: condition !! Assertion value.
      character(len=*), intent(in) :: label !! Human-readable assertion label.
      if (.not. condition) then
         print *, "FAIL: ", trim(label)
         error stop 1
      end if
   end subroutine check

   subroutine check_close(value, target, tolerance, label)
      !! Stops the test program when a scalar differs from its target beyond tolerance.
      real(dp), intent(in) :: value !! Computed value.
      real(dp), intent(in) :: target !! Expected value.
      real(dp), intent(in) :: tolerance !! Maximum permitted absolute error.
      character(len=*), intent(in) :: label !! Human-readable assertion label.
      if (abs(value - target) > tolerance) then
         print *, "FAIL: ", trim(label), value, target
         error stop 1
      end if
   end subroutine check_close

end program test_irregular
