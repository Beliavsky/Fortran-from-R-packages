program test_imputets
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
   use imputets_api, only : dp, na_interpolation, na_kalman, na_locf, na_ma, na_mean, na_random, &
                           na_remove, na_replace, na_seadec, na_seasplit, na_stats_result, stats_na
   implicit none

   integer :: failures

   failures = 0
   call test_replace_remove(failures)
   call test_mean(failures)
   call test_locf(failures)
   call test_ma(failures)
   call test_interpolation(failures)
   call test_random(failures)
   call test_kalman(failures)
   call test_seasonal(failures)
   call test_stats(failures)
   call test_maxgap(failures)

   if (failures /= 0) then
      write (*, '(a,i0)') 'FAILED tests: ', failures
      error stop 1
   end if
   write (*, '(a)') 'All imputeTS tests passed.'

contains

   subroutine test_replace_remove(failures)
      integer, intent(inout) :: failures !! Running number of failed assertions, incremented in place.
      real(dp), allocatable :: y(:)
      real(dp) :: x(5)

      x = [1.0_dp, nan(), 3.0_dp, nan(), 5.0_dp]
      y = na_replace(x, 9.0_dp)
      call check_close(y, [1.0_dp, 9.0_dp, 3.0_dp, 9.0_dp, 5.0_dp], 0.0_dp, 'na_replace', failures)
      y = na_remove(x)
      call check_close(y, [1.0_dp, 3.0_dp, 5.0_dp], 0.0_dp, 'na_remove', failures)
   end subroutine test_replace_remove

   subroutine test_mean(failures)
      integer, intent(inout) :: failures !! Running number of failed assertions, incremented in place.
      real(dp), allocatable :: y(:)
      real(dp) :: x(6)

      x = [1.0_dp, 2.0_dp, nan(), 2.0_dp, 5.0_dp, nan()]
      y = na_mean(x)
      call check_close(y, [1.0_dp, 2.0_dp, 2.5_dp, 2.0_dp, 5.0_dp, 2.5_dp], 1.0e-13_dp, 'mean', failures)
      y = na_mean(x, 'median')
      call check_scalar(y(3), 2.0_dp, 1.0e-13_dp, 'median', failures)
      y = na_mean(x, 'mode')
      call check_scalar(y(3), 2.0_dp, 1.0e-13_dp, 'mode', failures)

      x = [1.0_dp, 4.0_dp, nan(), 16.0_dp, 4.0_dp, nan()]
      y = na_mean(x, 'geometric')
      call check_scalar(y(3), 4.0_dp, 1.0e-12_dp, 'geometric mean', failures)
      y = na_mean(x, 'harmonic')
      call check_scalar(y(3), 2.56_dp, 1.0e-12_dp, 'harmonic mean', failures)
   end subroutine test_mean

   subroutine test_locf(failures)
      integer, intent(inout) :: failures !! Running number of failed assertions, incremented in place.
      real(dp), allocatable :: y(:)
      real(dp) :: x(7)

      x = [nan(), 2.0_dp, nan(), nan(), 5.0_dp, nan(), nan()]
      y = na_locf(x)
      call check_close(y, [2.0_dp, 2.0_dp, 2.0_dp, 2.0_dp, 5.0_dp, 5.0_dp, 5.0_dp], &
                       0.0_dp, 'locf reverse edge fill', failures)
      y = na_locf(x, 'nocb', 'keep')
      call check_true(ieee_is_nan(y(6)) .and. ieee_is_nan(y(7)), 'nocb keep trailing NaNs', failures)
      call check_scalar(y(1), 2.0_dp, 0.0_dp, 'nocb leading', failures)
   end subroutine test_locf

   subroutine test_ma(failures)
      integer, intent(inout) :: failures !! Running number of failed assertions, incremented in place.
      real(dp), allocatable :: y(:)
      real(dp) :: x(5)

      x = [1.0_dp, 2.0_dp, nan(), 4.0_dp, 5.0_dp]
      y = na_ma(x, 1, 'simple')
      call check_scalar(y(3), 3.0_dp, 1.0e-13_dp, 'moving average simple', failures)
      y = na_ma(x, 2, 'linear')
      call check_scalar(y(3), 3.0_dp, 1.0e-13_dp, 'moving average linear', failures)
      y = na_ma(x, 2, 'exponential')
      call check_scalar(y(3), 3.0_dp, 1.0e-13_dp, 'moving average exponential', failures)
   end subroutine test_ma

   subroutine test_interpolation(failures)
      integer, intent(inout) :: failures !! Running number of failed assertions, incremented in place.
      real(dp), allocatable :: y(:)
      real(dp) :: x(7)

      x = [nan(), 2.0_dp, nan(), 4.0_dp, nan(), 6.0_dp, nan()]
      y = na_interpolation(x, 'linear')
      call check_close(y, [2.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, 6.0_dp, 6.0_dp], &
                       1.0e-13_dp, 'linear interpolation', failures)
      y = na_interpolation(x, 'stine')
      call check_true(.not. any_nan(y), 'stine fills all gaps including edges', failures)
      call check_scalar(y(4), 4.0_dp, 0.0_dp, 'stine preserves observed values', failures)
      y = na_interpolation(x, 'spline')
      call check_true(.not. any_nan(y), 'spline fills all gaps', failures)
      call check_scalar(y(6), 6.0_dp, 0.0_dp, 'spline preserves observed values', failures)
   end subroutine test_interpolation

   subroutine test_random(failures)
      integer, intent(inout) :: failures !! Running number of failed assertions, incremented in place.
      integer, allocatable :: seed(:)
      integer :: seed_size
      real(dp), allocatable :: y(:)
      real(dp) :: x(5)

      call random_seed(size = seed_size)
      allocate(seed(seed_size))
      seed = 314159
      call random_seed(put = seed)
      x = [1.0_dp, nan(), 2.0_dp, nan(), 3.0_dp]
      y = na_random(x, -2.0_dp, 4.0_dp)
      call check_true(y(2) >= -2.0_dp .and. y(2) < 4.0_dp, 'random first bound', failures)
      call check_true(y(4) >= -2.0_dp .and. y(4) < 4.0_dp, 'random second bound', failures)
      call check_scalar(y(3), 2.0_dp, 0.0_dp, 'random preserves observation', failures)
   end subroutine test_random

   subroutine test_kalman(failures)
      integer, intent(inout) :: failures !! Running number of failed assertions, incremented in place.
      real(dp), allocatable :: y(:)
      real(dp) :: x(10)

      x = [1.0_dp, 2.0_dp, nan(), 4.0_dp, 5.0_dp, nan(), 7.0_dp, 8.0_dp, nan(), 10.0_dp]
      y = na_kalman(x, 'StructTS', .true.)
      call check_true(.not. any_nan(y), 'StructTS smoother fills gaps', failures)
      call check_scalar(y(5), 5.0_dp, 0.0_dp, 'StructTS preserves observations', failures)
      y = na_kalman(x, 'auto.arima', .false.)
      call check_true(.not. any_nan(y), 'auto.arima filtered fills gaps', failures)
   end subroutine test_kalman

   subroutine test_seasonal(failures)
      integer, intent(inout) :: failures !! Running number of failed assertions, incremented in place.
      real(dp), allocatable :: y(:)
      real(dp) :: x(16)

      x = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 1.1_dp, nan(), 3.1_dp, 4.1_dp, &
           1.2_dp, 2.2_dp, nan(), 4.2_dp, 1.3_dp, 2.3_dp, 3.3_dp, 4.3_dp]
      y = na_seasplit(x, 'interpolation', period = 4)
      call check_true(.not. any_nan(y), 'seasonal split fills gaps', failures)
      call check_scalar(y(1), 1.0_dp, 0.0_dp, 'seasonal split preserves observations', failures)
      y = na_seadec(x, 'interpolation', period = 4)
      call check_true(.not. any_nan(y), 'seasonal decomposition fills gaps', failures)
      call check_scalar(y(16), 4.3_dp, 0.0_dp, 'seasonal decomposition preserves observations', failures)

      x = [1.0_dp, nan(), nan(), 4.0_dp, 5.0_dp, 6.0_dp, 7.0_dp, 8.0_dp, &
           9.0_dp, 10.0_dp, 11.0_dp, 12.0_dp, 13.0_dp, 14.0_dp, 15.0_dp, 16.0_dp]
      y = na_seadec(x, 'interpolation', period = 1, maxgap = 0)
      call check_true(.not. any_nan(y), 'seasonal fallback preserves upstream early-return maxgap behavior', failures)
   end subroutine test_seasonal

   subroutine test_stats(failures)
      integer, intent(inout) :: failures !! Running number of failed assertions, incremented in place.
      type(na_stats_result) :: s
      real(dp) :: x(10)

      x = [1.0_dp, nan(), nan(), 2.0_dp, nan(), 3.0_dp, nan(), nan(), nan(), 4.0_dp]
      s = stats_na(x)
      call check_int(s%length_series, 10, 'stats length', failures)
      call check_int(s%number_nas, 6, 'stats NA count', failures)
      call check_int(s%number_na_gaps, 3, 'stats gap count', failures)
      call check_scalar(s%average_size_na_gaps, 2.0_dp, 1.0e-13_dp, 'stats average gap', failures)
      call check_scalar(s%percentage_nas, 60.0_dp, 1.0e-13_dp, 'stats percent', failures)
      call check_int(s%longest_na_gap, 3, 'stats longest gap', failures)
      call check_int(s%most_frequent_na_gap, 3, 'stats tie chooses largest', failures)
      call check_int(s%most_weighty_na_gap, 3, 'stats weightiest gap', failures)
      call check_int(s%distribution_na_gaps(1), 1, 'stats gaps of one', failures)
      call check_int(s%distribution_na_gaps(2), 1, 'stats gaps of two', failures)
      call check_int(s%distribution_na_gaps(3), 1, 'stats gaps of three', failures)
   end subroutine test_stats

   subroutine test_maxgap(failures)
      integer, intent(inout) :: failures !! Running number of failed assertions, incremented in place.
      real(dp), allocatable :: y(:)
      real(dp) :: x(8)

      x = [1.0_dp, nan(), 3.0_dp, nan(), nan(), nan(), 7.0_dp, 8.0_dp]
      y = na_interpolation(x, maxgap = 1)
      call check_scalar(y(2), 2.0_dp, 1.0e-13_dp, 'maxgap keeps short imputation', failures)
      call check_true(all_nan(y(4:6)), 'maxgap restores long gap', failures)
   end subroutine test_maxgap


   real(dp) function nan() result(value)
      use, intrinsic :: ieee_arithmetic, only : ieee_quiet_nan, ieee_value
      value = ieee_value(0.0_dp, ieee_quiet_nan)
   end function nan

   logical function any_nan(x) result(found)
      real(dp), intent(in) :: x(:) !! Numeric series checked for at least one NaN.
      integer :: i

      found = .false.
      do i = 1, size(x)
         if (ieee_is_nan(x(i))) then
            found = .true.
            return
         end if
      end do
   end function any_nan

   logical function all_nan(x) result(found)
      real(dp), intent(in) :: x(:) !! Numeric series checked to determine whether every element is NaN.
      integer :: i

      found = .true.
      do i = 1, size(x)
         if (.not. ieee_is_nan(x(i))) then
            found = .false.
            return
         end if
      end do
   end function all_nan

   subroutine check_close(actual, expected, tolerance, label, failures)
      real(dp), intent(in) :: actual(:) !! Computed vector being checked.
      real(dp), intent(in) :: expected(:) !! Expected vector with the same shape as actual.
      real(dp), intent(in) :: tolerance !! Maximum allowed absolute elementwise error.
      character(len=*), intent(in) :: label !! Human-readable assertion label printed when the check fails.
      integer, intent(inout) :: failures !! Running number of failed assertions, incremented on failure.

      if (size(actual) /= size(expected)) then
         failures = failures + 1
         write (*, '(a)') 'FAIL size: ' // trim(label)
      else if (any(abs(actual - expected) > tolerance)) then
         failures = failures + 1
         write (*, '(a)') 'FAIL values: ' // trim(label)
      end if
   end subroutine check_close

   subroutine check_scalar(actual, expected, tolerance, label, failures)
      real(dp), intent(in) :: actual !! Computed scalar being checked.
      real(dp), intent(in) :: expected !! Expected scalar value.
      real(dp), intent(in) :: tolerance !! Maximum allowed absolute error.
      character(len=*), intent(in) :: label !! Human-readable assertion label printed when the check fails.
      integer, intent(inout) :: failures !! Running number of failed assertions, incremented on failure.

      if (abs(actual - expected) > tolerance) then
         failures = failures + 1
         write (*, '(a,2es24.14)') 'FAIL scalar: ' // trim(label) // ' ', actual, expected
      end if
   end subroutine check_scalar

   subroutine check_int(actual, expected, label, failures)
      integer, intent(in) :: actual !! Computed integer being checked.
      integer, intent(in) :: expected !! Expected integer value.
      character(len=*), intent(in) :: label !! Human-readable assertion label printed when the check fails.
      integer, intent(inout) :: failures !! Running number of failed assertions, incremented on failure.

      if (actual /= expected) then
         failures = failures + 1
         write (*, '(a,2i0)') 'FAIL integer: ' // trim(label) // ' ', actual, expected
      end if
   end subroutine check_int

   subroutine check_true(condition, label, failures)
      logical, intent(in) :: condition !! Assertion condition expected to be true.
      character(len=*), intent(in) :: label !! Human-readable assertion label printed when the check fails.
      integer, intent(inout) :: failures !! Running number of failed assertions, incremented on failure.

      if (.not. condition) then
         failures = failures + 1
         write (*, '(a)') 'FAIL condition: ' // trim(label)
      end if
   end subroutine check_true

end program test_imputets
