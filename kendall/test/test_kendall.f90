! SPDX-License-Identifier: GPL-2.0-or-later
program test_kendall
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, ieee_positive_inf, ieee_is_nan
   use kendall_module
   implicit none

   call test_known_examples()
   call test_exact_enumeration()
   call test_tied_variance()
   call test_seasons()
   call test_edge_cases()
   print '(a)', 'All Kendall tests passed.'

contains

   subroutine test_known_examples()
      real(dp) :: x(12), y(12), a(17), b(17), expected_variance
      type(kendall_result) :: result, reverse, swapped

      x = [1.5_dp, 1.5_dp, 3.0_dp, 4.0_dp, 6.0_dp, 6.0_dp, 6.0_dp, 8.0_dp, 9.5_dp, 9.5_dp, 11.0_dp, 12.0_dp]
      y = [2.5_dp, 2.5_dp, 7.0_dp, 4.5_dp, 1.0_dp, 4.5_dp, 6.0_dp, 11.5_dp, 11.5_dp, 8.5_dp, 8.5_dp, 10.0_dp]
      call kendall(x, y, result)
      call require(result%status == kendall_ok, 'Kendall book example succeeds')
      call close(result%s, 34.0_dp, 0.0_dp, 'book score')
      call close(result%d, sqrt(61.0_dp * 62.0_dp), 1.0e-13_dp, 'book denominator')
      expected_variance = 6709.0_dp / 33.0_dp
      call close(result%vars, expected_variance, 1.0e-12_dp, 'book variance')
      call close(result%sl, erfc(33.0_dp / sqrt(2.0_dp * expected_variance)), 1.0e-15_dp, 'book p-value')
      call kendall(x, -y, reverse)
      call close(reverse%tau, -result%tau, 1.0e-15_dp, 'sign reversal')
      call close(reverse%sl, result%sl, 1.0e-15_dp, 'two-sided symmetry')
      call kendall(y, x, swapped)
      call close(swapped%tau, result%tau, 1.0e-15_dp, 'variable symmetry')
      call close(swapped%vars, result%vars, 1.0e-13_dp, 'symmetric tie variance')
      a = [2.5_dp, 2.5_dp, 2.5_dp, 2.5_dp, 5.0_dp, 6.5_dp, 6.5_dp, &
         10.0_dp, 10.0_dp, 10.0_dp, 10.0_dp, 10.0_dp, 14.0_dp, 14.0_dp, 14.0_dp, 16.0_dp, 17.0_dp]
      b = [1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 2.0_dp, 1.0_dp, 1.0_dp, &
         2.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 2.0_dp, 2.0_dp, 2.0_dp]
      call kendall(a, b, result)
      call close(result%s, 34.0_dp, 0.0_dp, 'second book score')
      call close(result%d, sqrt(116.0_dp * 60.0_dp), 1.0e-13_dp, 'second book denominator')
   end subroutine test_known_examples

   subroutine test_exact_enumeration()
      integer :: n, total_pairs, k, j, total, score, p(8), counts(0:28), representatives(8,0:28)
      real(dp) :: expected, x(8), y(8)
      logical :: more
      type(kendall_result) :: result

      do n = 3, 8
         counts = 0
         representatives = 0
         p(1:n) = [(j, j = 1, n)]
         total_pairs = n * (n - 1) / 2
         total = 0
         do
            k = inversions(p(1:n))
            counts(k) = counts(k) + 1
            representatives(1:n,k) = p(1:n)
            total = total + 1
            call next_permutation(p(1:n), more)
            if (.not. more) exit
         end do
         x(1:n) = [(real(j, dp), j = 1, n)]
         do k = 0, total_pairs
            y(1:n) = real(representatives(1:n,k), dp)
            score = total_pairs - 2 * k
            expected = 0.0_dp
            do j = 0, total_pairs
               if (abs(total_pairs - 2 * j) >= abs(score)) expected = expected + real(counts(j), dp)
            end do
            expected = expected / real(total, dp)
            call kendall(x(1:n), y(1:n), result, kendall_exact)
            call require(result%status == kendall_ok .and. result%exact, 'exact enumeration status')
            call close(result%s, real(score, dp), 0.0_dp, 'enumerated concordance score')
            call close(result%sl, expected, 1.0e-14_dp, 'exhaustive permutation tail')
         end do
      end do
      x(1:8) = [(real(j, dp), j = 1, 8)]
      call kendall(x, x, result, kendall_exact)
      call close(result%sl, 2.0_dp / 40320.0_dp, 1.0e-18_dp, 'inclusive exact extreme tail')
   end subroutine test_exact_enumeration

   subroutine test_tied_variance()
      real(dp) :: x(7), y(7), permuted(7), sum_s, sum_s2, score
      integer :: p(7), i, j, case_id, total, non_tied_x, non_tied_y
      logical :: more
      type(kendall_result) :: result

      do case_id = 1, 3
         x = [0.0_dp, 0.0_dp, 0.0_dp, 2.0_dp, 2.0_dp, 4.0_dp, 5.0_dp]
         y = [1.0_dp, 1.0_dp, 3.0_dp, 3.0_dp, 3.0_dp, 5.0_dp, 5.0_dp]
         if (case_id == 2) x = [(real(i, dp), i = 1, 7)]
         if (case_id == 3) y = [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp]
         p = [(i, i = 1, 7)]
         sum_s = 0.0_dp
         sum_s2 = 0.0_dp
         total = 0
         do
            permuted = y(p)
            score = 0.0_dp
            do i = 1, 6
               do j = i + 1, 7
                  score = score + signum(x(j) - x(i)) * signum(permuted(j) - permuted(i))
               end do
            end do
            sum_s = sum_s + score
            sum_s2 = sum_s2 + score * score
            total = total + 1
            call next_permutation(p, more)
            if (.not. more) exit
         end do
         call kendall(x, y, result)
         call close(sum_s / real(total, dp), 0.0_dp, 0.0_dp, 'permutation null mean')
         call close(result%vars, sum_s2 / real(total, dp), 1.0e-13_dp, 'tie variance by enumeration')
         non_tied_x = 0
         non_tied_y = 0
         do i = 1, 6
            do j = i + 1, 7
               if (x(i) /= x(j)) non_tied_x = non_tied_x + 1
               if (y(i) /= y(j)) non_tied_y = non_tied_y + 1
            end do
         end do
         call close(result%d, sqrt(real(non_tied_x * non_tied_y, dp)), 1.0e-13_dp, 'tie denominator enumeration')
      end do
   end subroutine test_tied_variance

   subroutine test_seasons()
      real(dp) :: x(18), score, denom, variance, missing(18), nan
      integer :: i, season
      type(kendall_result) :: result, part, with_missing, shifted

      x = [1.0_dp, 5.0_dp, 8.0_dp, 2.0_dp, 5.0_dp, 6.0_dp, 3.0_dp, 5.0_dp, &
         7.0_dp, 4.0_dp, 5.0_dp, 3.0_dp, 5.0_dp, 5.0_dp, 4.0_dp, 6.0_dp, 5.0_dp, 2.0_dp]
      score = 0.0_dp
      denom = 0.0_dp
      variance = 0.0_dp
      do season = 1, 3
         call mann_kendall(x(season::3), part)
         score = score + part%s
         denom = denom + part%d
         variance = variance + part%vars
      end do
      call seasonal_mann_kendall(x, 3, result)
      call require(result%status == kendall_ok, 'constant seasonal component accepted')
      call close(result%s, score, 0.0_dp, 'seasonal score aggregation')
      call close(result%d, denom, 0.0_dp, 'seasonal denominator aggregation')
      call close(result%vars, variance, 0.0_dp, 'seasonal variance aggregation')
      call close(result%sl, erfc(abs(score) / sqrt(2.0_dp * variance)), 1.0e-15_dp, 'seasonal uncorrected normal')
      call seasonal_mann_kendall([(real(i, dp), i = 1, 12)], 3, shifted)
      call close(shifted%s, 18.0_dp, 0.0_dp, 'monotonic seasonal score')
      call close(shifted%vars, 26.0_dp, 1.0e-14_dp, 'monotonic seasonal variance')
      call close(shifted%tau, 1.0_dp, 1.0e-15_dp, 'monotonic seasonal tau')
      nan = ieee_value(0.0_dp, ieee_quiet_nan)
      missing = x
      missing(4) = nan
      call seasonal_mann_kendall(missing, 3, with_missing)
      call require(with_missing%status == kendall_ok .and. with_missing%n_used == 17, 'seasonal slot-preserving missingness')
      call seasonal_mann_kendall(x(1:10), 3, shifted)
      call require(shifted%status == kendall_ok, 'partial last cycle is allowed')
      call seasonal_mann_kendall(missing(1:10), 3, shifted)
      call require(shifted%status == kendall_ok, 'unequal season lengths after omission')
      missing(7) = nan
      call seasonal_mann_kendall(missing(1:10), 3, shifted)
      call require(shifted%status == kendall_too_few, 'short complete season rejected')
      call seasonal_mann_kendall(x, 0, shifted)
      call require(shifted%status == kendall_invalid_period, 'zero period rejected')
      call seasonal_mann_kendall(x, -1, shifted)
      call require(shifted%status == kendall_invalid_period, 'negative period rejected')
      call seasonal_mann_kendall(x, 100, shifted)
      call require(shifted%status == kendall_too_few, 'oversized period rejected')
      call seasonal_mann_kendall([1.0_dp, 2.0_dp, 1.0_dp, 2.0_dp, 1.0_dp, 2.0_dp], 2, shifted)
      call require(shifted%status == kendall_degenerate .and. ieee_is_nan(shifted%tau), 'all constant seasons')
   end subroutine test_seasons

   subroutine test_edge_cases()
      real(dp) :: nan, inf, large, tiny_offset(5), x(101)
      type(kendall_result) :: result, reference
      integer :: i

      nan = ieee_value(0.0_dp, ieee_quiet_nan)
      inf = ieee_value(0.0_dp, ieee_positive_inf)
      call kendall([1.0_dp, 2.0_dp], [1.0_dp], result)
      call require(result%status == kendall_size_mismatch, 'unequal lengths rejected')
      call mann_kendall([1.0_dp, 2.0_dp], result)
      call require(result%status == kendall_too_few .and. ieee_is_nan(result%s), 'short input rejected')
      call mann_kendall([real(dp) ::], result)
      call require(result%status == kendall_too_few, 'empty input rejected')
      call mann_kendall([nan, nan, nan, nan], result)
      call require(result%status == kendall_too_few .and. result%n_used == 0, 'all missing input')
      call mann_kendall([1.0_dp, 2.0_dp, nan, nan], result)
      call require(result%status == kendall_too_few .and. result%n_used == 2, 'too few after omission')
      call kendall([1.0_dp, 2.0_dp, 3.0_dp, inf], [1.0_dp, 2.0_dp, 3.0_dp, nan], result)
      call require(result%status == kendall_ok .and. result%n_used == 3, 'NaN pair removed before finite check')
      call mann_kendall([1.0_dp, 2.0_dp, inf], result)
      call require(result%status == kendall_nonfinite, 'infinity rejected')
      call mann_kendall([1.0_dp, 1.0_dp, 1.0_dp], result)
      call require(result%status == kendall_degenerate .and. ieee_is_nan(result%sl), 'constant series undefined')
      call close(result%s + result%d + result%vars, 0.0_dp, 0.0_dp, 'constant contributions vanish')
      call mann_kendall([1.0_dp, 2.0_dp, 3.0_dp], result)
      call require(result%status == kendall_ok, 'three observations supported')
      call close(result%sl, erfc(2.0_dp / sqrt(22.0_dp / 3.0_dp)), 1.0e-15_dp, 'three-point normal p-value')
      call mann_kendall([1.0_dp, 2.0_dp, 3.0_dp], result, kendall_exact)
      call close(result%sl, 1.0_dp / 3.0_dp, 1.0e-15_dp, 'three-point exact p-value')
      call mann_kendall([1.0_dp, 1.0_dp, 3.0_dp], result, kendall_exact)
      call require(result%status == kendall_exact_ties .and. ieee_is_nan(result%sl), 'exact ties rejected')
      call mann_kendall([1.0_dp, 2.0_dp, 3.0_dp], result, -10)
      call require(result%status == kendall_invalid_method, 'unknown method rejected')
      x = [(real(i, dp), i = 1, 101)]
      call mann_kendall(x, result, kendall_exact)
      call require(result%status == kendall_exact_limit, 'exact resource limit explicit')
      call mann_kendall(x(1:100), result, kendall_exact)
      call close(result%sl / (2.0_dp / product(x(1:100))), 1.0_dp, 1.0e-13_dp, '100-point exact tail')
      call mann_kendall(x, result)
      call require(result%sl > 0.0_dp .and. result%sl < 1.0e-40_dp, 'small normal tail does not cancel to zero')
      call mann_kendall([4.0_dp, nan, 1.0_dp, 6.0_dp, 2.0_dp], result)
      call mann_kendall([4.0_dp, 1.0_dp, 6.0_dp, 2.0_dp], reference)
      call close(result%tau, reference%tau, 1.0e-15_dp, 'Mann-Kendall missingness')
      tiny_offset = 1.0_dp + [(real(i, dp) * 1.0e-10_dp, i = 1, 5)]
      call mann_kendall(tiny_offset, result)
      call close(result%tau, 1.0_dp, 1.0e-15_dp, 'dp distinctions preserved')
      large = huge(1.0_dp)
      call kendall([-large, 0.0_dp, large], [-large, 0.0_dp, large], result)
      call require(result%status == kendall_ok, 'score comparisons avoid subtraction overflow')
      call close(result%tau, 1.0_dp, 0.0_dp, 'extreme finite ordered values')
   end subroutine test_edge_cases

   pure integer function inversions(p) result(total)
      integer, intent(in) :: p(:) !! Permutation whose strictly inverted index pairs are counted independently.
      integer :: i, j

      total = 0
      do i = 1, size(p) - 1
         do j = i + 1, size(p)
            if (p(i) > p(j)) total = total + 1
         end do
      end do
   end function inversions

   pure elemental real(dp) function signum(value) result(signed)
      real(dp), intent(in) :: value !! Finite difference from small deterministic test data.

      signed = 0.0_dp
      if (value > 0.0_dp) signed = 1.0_dp
      if (value < 0.0_dp) signed = -1.0_dp
   end function signum

   pure subroutine next_permutation(p, more)
      integer, intent(inout) :: p(:) !! Ascending-lexicographic permutation, advanced in-place when possible.
      logical, intent(out) :: more !! True if a successor exists; false at the last permutation.
      integer :: i, j, swap, first, last

      i = size(p) - 1
      do while (i >= 1)
         if (p(i) < p(i + 1)) exit
         i = i - 1
      end do
      more = i >= 1
      if (.not. more) return
      j = size(p)
      do while (p(j) <= p(i))
         j = j - 1
      end do
      swap = p(i)
      p(i) = p(j)
      p(j) = swap
      first = i + 1
      last = size(p)
      do while (first < last)
         swap = p(first)
         p(first) = p(last)
         p(last) = swap
         first = first + 1
         last = last - 1
      end do
   end subroutine next_permutation

   subroutine close(actual, expected, tolerance, description)
      real(dp), intent(in) :: actual !! Computed statistic, required finite for this assertion.
      real(dp), intent(in) :: expected !! Independently specified finite reference statistic.
      real(dp), intent(in) :: tolerance !! Nonnegative absolute error allowed by the assertion.
      character(len = *), intent(in) :: description !! Diagnostic text identifying the tested invariant.

      if (ieee_is_nan(actual) .or. abs(actual - expected) > tolerance) then
         print *, description, ': actual = ', actual, ', expected = ', expected
         error stop 1
      end if
   end subroutine close

   subroutine require(condition, description)
      logical, intent(in) :: condition !! Logical invariant that must hold for this deterministic test.
      character(len = *), intent(in) :: description !! Diagnostic explaining the invariant if it fails.

      if (.not. condition) then
         print *, description
         error stop 1
      end if
   end subroutine require

end program test_kendall
