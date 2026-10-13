! SPDX-License-Identifier: GPL-3.0-or-later
program test_r_reference
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use infotheo
   use r_reference_data
   implicit none
   integer :: i, failures, status
   integer, allocatable :: bins(:, :)
   character(len=6), parameter :: methods(4) = [character(len=6) :: 'emp', 'mm', 'sg_r', 'shrink']
   failures = 0
   do i = 1, 4
      call check_measures(reference_data, trim(methods(i)), reference_values(:, i))
      call check_measures(missing_data, trim(methods(i)), missing_values(:, i))
   end do
   call discretize(continuous_data, bins, 'equalfreq', 5, status)
   call check(status == infotheo_ok .and. all(bins == ef_bins), 'R equal-frequency fixture')
   call discretize(continuous_data, bins, 'equalwidth', 5, status)
   call check(status == infotheo_ok .and. all(bins == ew_bins), 'R equal-width fixture')
   call discretize(continuous_data, bins, 'globalequalwidth', 5, status)
   call check(status == infotheo_ok .and. all(bins == global_bins), 'R global-width fixture')
   if (failures /= 0) error stop 'R reference regression failure'
   print '(a)', 'All 144 actual-R information values and 72 discretization labels passed.'
contains
   subroutine check_measures(x, method, reference)
      integer, intent(in) :: x(:, :) !! Twelve observations by three discrete variables, possibly missing.
      character(len=*), intent(in) :: method !! Fortran estimator corresponding to the R fixture's method.
      real(dp), intent(in) :: reference(:) !! Eighteen actual-R values in the generator's documented order.
      real(dp) :: actual(18)
      real(dp), allocatable :: mi(:, :)
      integer :: j
      actual(1) = entropy(x, method)
      do j = 1, 3
         actual(j + 1) = entropy(x(:, j), method)
      end do
      actual(5) = condentropy(x(:, 1:2), x(:, 3:3), method)
      actual(6) = mutinformation(x(:, 1:1), x(:, 2:3), method)
      actual(7) = condinformation(x(:, 1:1), x(:, 2:2), x(:, 3:3), method)
      actual(8) = multiinformation(x, method)
      actual(9) = interinformation(x, method)
      mi = mutinformation_matrix(x, method)
      actual(10:18) = reshape(mi, [9])
      do j = 1, 18
         call check(ieee_is_finite(actual(j)), 'finite R-reference result ' // method)
         call check(abs(actual(j) - reference(j)) < 5.0e-12_dp, 'R information fixture ' // method)
         if (abs(actual(j) - reference(j)) >= 5.0e-12_dp) then
            print '(a,i0,2es25.16)', 'Index / Fortran / R: ', j, actual(j), reference(j)
         end if
      end do
   end subroutine check_measures
   subroutine check(condition, label)
      logical, intent(in) :: condition !! True when the current R-reference comparison passes.
      character(len=*), intent(in) :: label !! Regression description printed if a comparison fails.
      if (.not. condition) then
         print '(a)', 'FAIL: ' // label
         failures = failures + 1
      end if
   end subroutine check
end program test_r_reference
