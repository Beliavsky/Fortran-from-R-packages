! SPDX-License-Identifier: GPL-3.0-or-later
program test_infotheo
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_value, ieee_quiet_nan, ieee_positive_inf
   use infotheo
   implicit none
   integer :: failures, status, i, j
   integer :: xor_data(4, 3), independent(4, 2), counts(10)
   integer :: with_missing(4, 2), same(1, 1), empty(0, 2), no_columns(4, 0), too_many(1, 21)
   integer, allocatable :: bins(:, :), copied(:, :), large(:, :)
   real(dp), allocatable :: matrix(:, :)
   real(dp) :: continuous(8, 2), bounds(5, 1), nan_data(4, 2), all_missing(2, 1)
   real(dp) :: huge_data(3, 1), degenerate(4, 1), h, log2
   character(len=6), parameter :: methods(5) = [character(len=6) :: 'emp', 'mm', 'sg', 'sg_r', 'shrink']
   real(dp), parameter :: counts_expected(5) = [0.8018185525433373_dp, 0.9018185525433373_dp, &
      0.7730536807232604_dp, 0.7139682539682545_dp, 0.8625257778161765_dp]

   failures = 0
   log2 = log(2.0_dp)
   xor_data(:, 1) = [0, 0, 1, 1]
   xor_data(:, 2) = [0, 1, 0, 1]
   xor_data(:, 3) = [0, 1, 1, 0]
   independent = xor_data(:, :2)
   counts = [0, 1, 1, 2, 2, 2, 2, 2, 2, 2]
   do i = 1, size(methods)
      call close(entropy(counts, trim(methods(i))), counts_expected(i), 3.0e-12_dp, 'nonuniform estimator')
      same = 123
      call close(entropy(same, trim(methods(i))), 0.0_dp, 0.0_dp, 'singleton estimator')
      call close(entropy([1, 1, 1], trim(methods(i))), 0.0_dp, 0.0_dp, 'constant estimator')
   end do
   call close(entropy(xor_data), 2.0_dp * log2, 2.0e-14_dp, 'XOR joint entropy')
   call close(multiinformation(xor_data), log2, 2.0e-14_dp, 'XOR total correlation')
   call close(interinformation(xor_data), -log2, 2.0e-14_dp, 'XOR co-information sign')
   call close(mutinformation(xor_data(:, 1:1), xor_data(:, 2:2)), 0.0_dp, 1.0e-14_dp, 'independence')
   call close(condinformation(xor_data(:, 1:1), xor_data(:, 2:2), xor_data(:, 3:3)), &
      log2, 2.0e-14_dp, 'conditional XOR dependence')
   call close(condinformation(xor_data(:, 1:1), xor_data(:, 2:2), xor_data(:, 2:3)), &
      0.0_dp, 2.0e-14_dp, 'multicolumn conditioning includes Y')
   call close(condinformation(xor_data(:, :2), xor_data(:, 3:3), xor_data(:, 1:1)), &
      log2, 2.0e-14_dp, 'multicolumn target block')
   call close(condentropy(xor_data(:, :2), xor_data(:, 3:3)), log2, 2.0e-14_dp, 'block conditional entropy')
   call close(condentropy(xor_data), entropy(xor_data), 0.0_dp, 'absent conditioning')
   call close(condinformation(xor_data(:, :2), xor_data(:, 3:3)), &
      mutinformation(xor_data(:, :2), xor_data(:, 3:3)), 0.0_dp, 'absent conditional MI')
   call close(interinformation(xor_data(:, 1:1)), log2, 1.0e-14_dp, 'one-variable interaction')
   call close(multiinformation(xor_data(:, 1:1)), 0.0_dp, 1.0e-14_dp, 'one-variable total correlation')
   call close(interinformation(independent), 0.0_dp, 1.0e-14_dp, 'two-variable interaction')
   matrix = mutinformation_matrix(xor_data)
   do i = 1, 3
      call close(matrix(i, i), log2, 1.0e-14_dp, 'MI diagonal entropy')
      do j = 1, i - 1
         call close(matrix(i, j), 0.0_dp, 1.0e-14_dp, 'MI off-diagonal independence')
         call close(matrix(i, j), matrix(j, i), 0.0_dp, 'MI symmetry')
      end do
   end do
   matrix = mutinformation_matrix(independent, 'mm')
   call close(matrix(1, 2), -0.125_dp, 1.0e-14_dp, 'MM matrix retains negative estimate')
   call close(mutinformation(independent(:, 1:1), independent(:, 2:2), 'mm'), &
      0.0_dp, 0.0_dp, 'MM scalar clips negative estimate')
   call close(natstobits(log2), 1.0_dp, 1.0e-15_dp, 'exact conversion')
   call check(all(abs(natstobits([0.0_dp, log2, 2.0_dp * log2]) - [0.0_dp, 1.0_dp, 2.0_dp]) &
      < 1.0e-14_dp), 'elemental conversion')

   with_missing = independent
   with_missing(2, 2) = na_integer
   with_missing(4, 1) = na_integer
   call close(entropy(with_missing), log2, 1.0e-14_dp, 'joint complete observations')
   call close(entropy(with_missing(:, 1)), entropy([0, 0, 1]), 1.0e-14_dp, 'marginal complete observations')
   call close(condentropy(with_missing(:, 1:1), with_missing(:, 2:2)), &
      entropy(with_missing) - entropy(with_missing(:, 2)), 1.0e-14_dp, 'available-case R convention')
   call check(ieee_is_nan(entropy([na_integer, na_integer])), 'all-missing entropy invalid')
   call check(ieee_is_nan(entropy(empty)), 'empty entropy invalid')
   call check(ieee_is_nan(entropy(counts, 'bad')), 'unknown estimator invalid')
   call check(ieee_is_nan(condentropy(independent, same)), 'conditional shape mismatch invalid')
   call check(ieee_is_nan(mutinformation(independent, same)), 'mutual shape mismatch invalid')
   call check(ieee_is_nan(condinformation(independent, independent, same)), 'conditioning shape mismatch invalid')
   too_many = 0
   call check(ieee_is_nan(interinformation(too_many)), 'exponential interaction guard')
   call close(entropy(no_columns), 0.0_dp, 0.0_dp, 'empty variable set entropy')
   call close(interinformation(no_columns), 0.0_dp, 0.0_dp, 'empty variable set interaction')
   call close(multiinformation(no_columns), 0.0_dp, 0.0_dp, 'empty variable set total correlation')
   allocate(large(100000, 1))
   large(:90000, 1) = 1
   large(90001:, 1) = 2
   h = entropy(large, 'shrink')
   call check(h > 0.32_dp .and. h < 0.34_dp, 'shrink avoids squared integer overflow')
   deallocate(large)

   continuous(:, 1) = [(real(i, dp), i = 1, 8)]
   continuous(:, 2) = continuous(8:1:-1, 1)
   call discretize(continuous, bins, nbins = 3, status = status)
   call check(status == infotheo_ok, 'equal-frequency status')
   call check(all(bins(:, 1) == [1, 1, 1, 2, 2, 3, 3, 3]), 'legacy equal-frequency remainder ranks')
   call check(all(bins(:, 2) == bins(8:1:-1, 1)), 'equal-frequency descending observations')
   continuous(:, 1) = [1.0_dp, 1.0_dp, 2.0_dp, 2.0_dp, 2.0_dp, 3.0_dp, 3.0_dp, 4.0_dp]
   call discretize(continuous, bins, nbins = 3)
   call check(all(bins(:, 1) == [1, 1, 1, 1, 1, 3, 3, 3]), 'equal-frequency ties stay together')
   bounds(:, 1) = [0.0_dp, 1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp]
   call discretize(bounds, bins, disc = 'equalwidth', nbins = 2, status = status)
   call check(status == infotheo_ok, 'equal-width status')
   call check(all(bins(:, 1) == [1, 1, 2, 2, 2]), 'equal-width left-closed boundaries')
   call discretize(bounds, bins, disc = 'globalequalwidth', nbins = 2, status = status)
   call check(status == infotheo_ok, 'global-width status')
   call check(all(bins(:, 1) == [1, 1, 1, 2, 2]), 'global-width right-closed boundaries')
   bounds = bounds + 1.0e12_dp
   call discretize(bounds, bins, disc = 'equalwidth', nbins = 2)
   call check(all(bins(:, 1) == [1, 1, 2, 2, 2]), 'equal-width extrema beyond int32')
   huge_data(:, 1) = [-huge(0.0_dp), 0.0_dp, huge(0.0_dp)]
   call discretize(huge_data, bins, disc = 'equalwidth', nbins = 2)
   call check(all(bins(:, 1) == [1, 2, 2]), 'equal-width huge finite range')
   call discretize(huge_data, bins, disc = 'globalequalwidth', nbins = 2)
   call check(all(bins(:, 1) == [1, 1, 2]), 'global-width huge finite range')
   call discretize(huge_data, bins, disc = 'equalwidth', nbins = huge(0))
   call check(bins(3, 1) == huge(0), 'equal-width max integer bins avoids overflow')
   degenerate = 7.0_dp
   call discretize(degenerate, bins, nbins = 3)
   call check(all(bins == 1), 'constant equal-frequency bin')
   call discretize(degenerate, bins, 'equalwidth', 3)
   call check(all(bins == 1), 'constant equal-width bin')
   call discretize(degenerate, bins, 'globalequalwidth', 4)
   call check(all(bins == 2), 'constant global even bins')
   call discretize(degenerate, bins, 'globalequalwidth', 3)
   call check(all(bins == 2), 'constant global odd bins')
   degenerate = 0.0_dp
   call discretize(degenerate, bins, 'globalequalwidth', 4)
   call check(all(bins == 2), 'zero constant global bins')

   nan_data(:, 1) = [1.0_dp, ieee_value(0.0_dp, ieee_quiet_nan), 2.0_dp, 3.0_dp]
   nan_data(:, 2) = ieee_value(0.0_dp, ieee_quiet_nan)
   do i = 1, 3
      select case (i)
      case (1)
         call discretize(nan_data, bins, 'equalfreq', 2, status)
      case (2)
         call discretize(nan_data, bins, 'equalwidth', 2, status)
      case (3)
         call discretize(nan_data, bins, 'globalequalwidth', 2, status)
      end select
      call check(status == infotheo_ok, 'missing discretization status')
      call check(bins(2, 1) == na_integer, 'missing observation sentinel')
      call check(all(bins(:, 2) == na_integer), 'all-missing column sentinel')
      call check(all(pack(bins(:, 1), bins(:, 1) /= na_integer) >= 1), 'nonmissing bins lower bound')
      call check(all(pack(bins(:, 1), bins(:, 1) /= na_integer) <= 2), 'nonmissing bins upper bound')
   end do
   all_missing = ieee_value(0.0_dp, ieee_quiet_nan)
   call discretize(all_missing, bins, 'globalequalwidth', 2, status)
   call check(status == infotheo_ok .and. all(bins == na_integer), 'all-missing global matrix')
   call discretize(nan_data, bins, 'equalfreq', 4, status)
   call check(status == infotheo_invalid .and. all(bins == na_integer), 'bin count exceeds valid observations')
   call discretize(bounds, bins, 'unknown', 2, status)
   call check(status == infotheo_invalid .and. all(bins == na_integer), 'unknown discretization method')
   call discretize(bounds, bins, nbins = 0, status = status)
   call check(status == infotheo_invalid, 'nonpositive bin count invalid')
   call discretize(bounds, bins, 'globalequalwidth', 1, status)
   call check(status == infotheo_invalid, 'single global bin invalid as in cut')
   nan_data(1, 1) = ieee_value(0.0_dp, ieee_positive_inf)
   call discretize(nan_data, bins, 'equalwidth', 2, status)
   call check(status == infotheo_invalid .and. all(bins == na_integer), 'infinite discretization input invalid')
   call discretize(continuous, bins)
   call discretize(continuous, copied, nbins = int(real(size(continuous, 1), dp) ** (1.0_dp / 3.0_dp)))
   call check(all(bins == copied), 'default bin count')

   if (failures /= 0) then
      print '(a,i0)', 'FAILURES: ', failures
      error stop 1
   end if
   print '(a)', 'All infotheo analytic, missing-data, and discretization tests passed.'
contains
   subroutine check(condition, label)
      logical, intent(in) :: condition !! Assertion that must hold for this deterministic regression.
      character(len=*), intent(in) :: label !! Human-readable identification printed on failure.
      if (.not. condition) then
         print '(a)', 'FAIL: ' // label
         failures = failures + 1
      end if
   end subroutine check
   subroutine close(actual, expected, tolerance, label)
      real(dp), intent(in) :: actual !! Computed finite result being checked.
      real(dp), intent(in) :: expected !! Independent analytic or independently evaluated reference.
      real(dp), intent(in) :: tolerance !! Allowed absolute error in the same units as the result.
      character(len=*), intent(in) :: label !! Human-readable test identification printed on failure.
      call check(.not. ieee_is_nan(actual), label // ' is finite')
      call check(abs(actual - expected) <= tolerance, label)
      if (abs(actual - expected) > tolerance) print '(2es25.16)', actual, expected
   end subroutine close
end program test_infotheo
