! SPDX-License-Identifier: GPL-2.0-or-later
program test_exact_rank
   use, intrinsic :: iso_fortran_env, only : int64
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_is_finite, ieee_value, ieee_quiet_nan
   use exact_rank_tests
   implicit none
   type(rank_test_result) :: result
   type(rank_distribution) :: distribution
   real(dp) :: x(5), y(4), a(6), expected3(3), expected2(2), expected5(5), expected7(7)
   real(dp) :: population(9), score, probability, total, nan, uniform(2, 12), sign_uniform(3, 8)
   real(dp), allocatable :: scores(:), draws(:)
   integer, allocatable :: ranks(:)
   character(9), parameter :: alternatives(3) = ['two.sided', 'less     ', 'greater  ']
   character(16), parameter :: kinds(7) = [character(16) :: 'Data', 'Wilcoxon', 'NormalQuantile', &
      'AnsariBradley', 'Median', 'Savage', 'ConSal']
   integer :: unit, status, ex, alt, i, j, k, m, n, mask, count_cases, checks, l, pair
   logical :: exact, paired
   checks = 0
   x = [1.0_dp, 2.0_dp, 2.0_dp, 3.0_dp, 5.0_dp]
   y = [0.0_dp, 2.0_dp, 4.0_dp, 4.0_dp]
   a = [-2.0_dp, -1.0_dp, 0.0_dp, 1.0_dp, 3.0_dp, 3.0_dp]
   open(newunit = unit, file = 'test/reference_values.dat', status = 'old', action = 'read')
   do ex = 1, 2
      exact = ex == 1
      do alt = 1, 3
         call wilcox_exact(x, result, y, alternative = trim(alternatives(alt)), exact = exact)
         call check_reference3()
         call wilcox_exact(a, result, alternative = trim(alternatives(alt)), exact = exact)
         call check_reference3()
         call ansari_exact(x, y, result, alternative = trim(alternatives(alt)), exact = exact)
         call check_reference2()
         call perm_test(x, result, y, alternative = trim(alternatives(alt)), exact = exact)
         call check_reference3()
         call perm_test(a, result, alternative = trim(alternatives(alt)), exact = exact)
         call check_reference3()
      end do
   end do
   do alt = 1, 3
      call wilcox_exact(x, result, y, alternative = trim(alternatives(alt)), exact = .true., conf_level = 0.8_dp)
      read(unit, *) expected3
      call assert(result%status == rank_ok .and. result%has_conf_int, 'two-sample interval status')
      call near(result%conf_int(1), expected3(1), 'two-sample interval lower')
      call near(result%conf_int(2), expected3(2), 'two-sample interval upper')
      call near(result%estimate, expected3(3), 'two-sample estimate')
      call wilcox_exact(a, result, alternative = trim(alternatives(alt)), exact = .true., conf_level = 0.8_dp)
      read(unit, *) expected3
      call assert(result%status == rank_ok .and. result%has_conf_int, 'signed interval status')
      call near(result%conf_int(1), expected3(1), 'signed interval lower')
      call near(result%conf_int(2), expected3(2), 'signed interval upper')
      call near(result%estimate, expected3(3), 'signed estimate')
   end do
   do ex = 1, 2
      do j = 1, 2
         if (j == 1) then
            call perm_test([-1.2_dp, 0.5_dp, 2.7_dp, 4.0_dp], result, [0.1_dp, 0.3_dp, 2.3_dp], exact = ex == 1)
         else
            call perm_test([-1.2_dp, 0.5_dp, 2.7_dp, 4.0_dp], result, [0.1_dp, 0.3_dp, 2.3_dp], &
               exact = ex == 1, tolerance = 0.01_dp)
         end if
         call check_reference2()
         if (j == 1) then
            call perm_test([-2.1_dp, -1.2_dp, 0.0_dp, 0.7_dp, 3.3_dp], result, exact = ex == 1)
         else
            call perm_test([-2.1_dp, -1.2_dp, 0.0_dp, 0.7_dp, 3.3_dp], result, exact = ex == 1, tolerance = 0.01_dp)
         end if
         call check_reference2()
      end do
   end do
   do ex = 1, 2
      do i = 1, 7
         call cscores([2.0_dp, 1.0_dp, 2.0_dp, 4.0_dp, 6.0_dp], scores, status, &
            trim(kinds(i)), integer_scores = ex == 2)
         call assert(status == rank_ok, 'score status')
         read(unit, *) expected5
         do j = 1, 5
            call near(scores(j), expected5(j), 'R scoring')
         end do
      end do
   end do
   do ex = 1, 2
      call cscores_surv([1.0_dp, 2.0_dp, 2.0_dp, 3.0_dp, 5.0_dp], [1, 0, 1, 1, 0], &
         scores, status, integer_scores = ex == 2)
      read(unit, *) expected5
      call assert(status == rank_ok, 'survival score status')
      do j = 1, 5
         call near(scores(j), expected5(j), 'R log-rank scoring')
      end do
   end do
   call permutation_distribution([-1.0_dp, 0.0_dp, 0.5_dp, 2.5_dp, 2.5_dp], 2, distribution, status)
   call assert(status == rank_ok, 'distribution status')
   read(unit, *) expected7
   call assert(size(distribution%support) == 7, 'distribution support size')
   do j = 1, 7
      call near(distribution%support(j), expected7(j), 'R distribution support')
   end do
   read(unit, *) expected7
   do j = 1, 7
      call near(dperm(distribution%support(j), distribution), expected7(j), 'R density')
   end do
   do alt = 1, 3
      read(unit, *) expected7
      do j = 1, 7
         select case (alt)
         case (1)
            probability = pperm(distribution%support(j), distribution, 'less')
         case (2)
            probability = pperm(distribution%support(j), distribution, 'greater')
         case (3)
            probability = pperm(distribution%support(j), distribution, 'two.sided')
         end select
         call near(probability, expected7(j), 'R inclusive tails')
      end do
   end do
   read(unit, *) expected7
   scores = qperm([0.0_dp, 0.1_dp, 0.25_dp, 0.5_dp, 0.75_dp, 0.99_dp, 1.0_dp], distribution)
   do j = 1, 7
      call near(scores(j), expected7(j), 'R quantiles')
   end do
   do alt = 1, 3
      call ansari_exact([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp], [1.0_dp, 2.0_dp, 4.0_dp, 6.0_dp, 9.0_dp], &
         result, trim(alternatives(alt)), exact = .true., conf_level = 0.8_dp)
      read(unit, *) expected3
      call assert(result%status == rank_ok .and. result%has_conf_int, 'Ansari interval status')
      call near(result%conf_int(1), expected3(1), 'Ansari interval lower')
      call near(result%conf_int(2), expected3(2), 'Ansari interval upper')
      call near(result%estimate, expected3(3), 'Ansari scale estimate')
   end do
   call wilcox_exact([-1.5_dp, -0.5_dp, -2.5_dp, 2.0_dp, -0.5_dp, 4.0_dp, 3.5_dp, 0.0_dp, &
      -3.0_dp, -0.5_dp, -3.5_dp], result, mu = 1.0_dp, exact = .true., conf_level = 0.9_dp)
   read(unit, *) expected3
   call assert(result%status == rank_ok .and. result%has_conf_int, 'shifted interval status')
   call near(result%conf_int(1), expected3(1), 'shifted interval lower')
   call near(result%conf_int(2), expected3(2), 'shifted interval upper')
   call near(result%estimate, expected3(3), 'shifted interval estimate')
   close(unit)
   ! Independent exhaustive enumeration, rather than a second DP implementation.
   do n = 1, 9
      population(:n) = [(real(modulo(7 * i + n, 11) - 5, dp) / 2.0_dp, i = 1, n)]
      do m = 1, n
         do pair = 0, 1
            if (pair == 1 .and. m /= n) cycle
            paired = pair == 1
            if (paired) then
               scores = abs(population(:n))
            else
               scores = population(:n)
            end if
            call permutation_distribution(scores, m, distribution, status, paired)
            call assert(status == rank_ok, 'enumerated distribution status')
            count_cases = 0
            do mask = 0, 2 ** n - 1
               if (.not. paired .and. popcnt(mask) /= m) cycle
               count_cases = count_cases + 1
            end do
            do j = 1, size(distribution%support)
               total = 0.0_dp
               do mask = 0, 2 ** n - 1
                  if (.not. paired .and. popcnt(mask) /= m) cycle
                  score = 0.0_dp
                  do i = 1, n
                     if (btest(mask, i - 1)) score = score + scores(i)
                  end do
                  if (score == distribution%support(j)) total = total + 1.0_dp
               end do
               call near(distribution%probability(j), total / real(count_cases, dp), 'exhaustive mass')
            end do
            call near(sum(distribution%probability), 1.0_dp, 'normalization')
         end do
      end do
   end do
   call permutation_distribution([7.0_dp, -1.5_dp, 9.0_dp, -3.5_dp, 1.0_dp, -6.0_dp, -2.5_dp], &
      3, distribution, status, factor = 2.0_dp)
   call near(pperm(2.0_dp, distribution, 'two.sided'), 34.0_dp / 35.0_dp, 'symmetric tail regression')
   call irank([3.0_dp, 1.0_dp, 3.0_dp, 2.0_dp], ranks, status)
   call assert(status == rank_ok .and. all(ranks == [4, 1, 4, 2]), 'integer maximum ranks')
   call irank([3.0_dp, 1.0_dp, 3.0_dp, 2.0_dp], ranks, status, [2, 4, 1, 3])
   call assert(status == rank_ok .and. all(ranks == [4, 1, 4, 2]), 'presorted integer ranks')
   call irank([1.0_dp, 2.0_dp], ranks, status, [1, 1])
   call assert(status == rank_invalid, 'duplicate order rejection')
   call cscores_factor([1, 2, 2, 1], scores, status)
   call assert(status == rank_ok .and. all(scores == [1.0_dp, 2.0_dp, 2.0_dp, 1.0_dp]), 'factor scores')
   call cscores_factor([1, 3], scores, status)
   call assert(status == rank_invalid, 'factor level rejection')
   k = 0
   do i = 0, 3
      do j = 0, 2
         k = k + 1
         uniform(:, k) = [(real(i, dp) + 0.25_dp) / 4.0_dp, (real(j, dp) + 0.25_dp) / 3.0_dp]
      end do
   end do
   call rperm([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp], 2, uniform, draws, status)
   call assert(status == rank_ok, 'sampling status')
   call near(sum(draws) / 12.0_dp, 5.0_dp, 'balanced sampling mean')
   call simulate_distribution([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp], 2, uniform, distribution, status)
   call assert(status == rank_ok, 'simulation status')
   call near(dperm(5.0_dp, distribution), 1.0_dp / 3.0_dp, 'balanced simulation mass')
   do j = 1, 8
      do i = 1, 3
         sign_uniform(i, j) = merge(0.25_dp, 0.75_dp, btest(j - 1, i - 1))
      end do
   end do
   call simulate_distribution([1.0_dp, 2.0_dp, 3.0_dp], 3, sign_uniform, distribution, status)
   call near(dperm(3.0_dp, distribution), 0.25_dp, 'sign-flip simulation')
   nan = ieee_value(0.0_dp, ieee_quiet_nan)
   call assert(ieee_is_nan(pperm(nan, distribution)), 'NaN CDF')
   call assert(ieee_is_nan(pperm(1.0_dp, distribution, 'bad')), 'bad alternative')
   call assert(ieee_is_nan(qperm(-0.1_dp, distribution)), 'bad quantile probability')
   call permutation_distribution([1.0_dp, nan], 1, distribution, status)
   call assert(status == rank_invalid, 'nonfinite input rejection')
   call permutation_distribution([1.0_dp, 2.0_dp], 1, distribution, status, max_cells = 1_int64)
   call assert(status == rank_limit, 'cell guard')
   call permutation_distribution([1.0e20_dp, 1.0_dp], 1, distribution, status)
   call assert(status == rank_limit, 'integer conversion guard')
   call permutation_distribution([1.0e308_dp, 1.0e308_dp, 0.0_dp], 1, distribution, status, factor = 1.0e-300_dp)
   call assert(status == rank_numeric, 'tiny factor expectation overflow guard')
   call permutation_distribution([1.0e308_dp, 1.0e308_dp], 2, distribution, status, &
      paired = .false., factor = 1.0e-300_dp)
   call assert(status == rank_numeric, 'tiny factor support overflow guard')
   call permutation_distribution([-1.0_dp, 1.0_dp], 2, distribution, status)
   call assert(status == rank_invalid, 'negative sign-flip rejection')
   call permutation_distribution([0.0_dp, 0.0_dp], 2, distribution, status)
   call assert(status == rank_ok, 'degenerate sign-flip status')
   call near(dperm(0.0_dp, distribution), 1.0_dp, 'degenerate sign flips')
   scores = [(1.0_dp, i = 1, 1100)]
   call permutation_distribution(scores, 1100, distribution, status)
   call assert(status == rank_ok, 'probability DP avoids count overflow')
   call near(sum(distribution%probability), 1.0_dp, 'large distribution normalized')
   call near(sum(distribution%support * distribution%probability), 550.0_dp, 'large distribution expectation')
   call wilcox_exact([1.0_dp, nan, 2.0_dp], result, [0.0_dp, 100.0_dp, 0.0_dp], paired = .true.)
   call assert(result%status == rank_ok, 'paired missing-value omission')
   call near(result%statistic, 3.0_dp, 'paired statistic')
   call wilcox_exact([1.0_dp], result, [1.0_dp, 2.0_dp], paired = .true.)
   call assert(result%status == rank_invalid, 'paired shape rejection')
   print '(a,i0,a)', 'exactRankTests: ', checks, ' deterministic checks passed.'
contains
   subroutine assert(condition, label)
      logical, intent(in) :: condition !! Test condition that must hold.
      character(*), intent(in) :: label !! Diagnostic identifying the failed assertion.
      checks = checks + 1
      if (.not. condition) then
         print *, 'FAIL: ', label
         error stop 1
      end if
   end subroutine assert

   subroutine near(actual, expected, label)
      real(dp), intent(in) :: actual !! Computed numerical result.
      real(dp), intent(in) :: expected !! Independent or upstream R reference value.
      character(*), intent(in) :: label !! Diagnostic identifying the numerical comparison.
      logical :: matched
      matched = .false.
      if (ieee_is_nan(expected)) then
         matched = ieee_is_nan(actual)
      else if (ieee_is_finite(expected)) then
         if (ieee_is_finite(actual)) matched = abs(actual - expected) <= 3.0e-12_dp * max(1.0_dp, abs(expected))
      else
         matched = actual == expected
      end if
      if (.not. matched) print *, label, ': actual=', actual, ' expected=', expected
      call assert(matched, label)
   end subroutine near

   subroutine check_reference3()
      read(unit, *) expected3
      call assert(result%status == rank_ok, 'R test status')
      call near(result%statistic, expected3(1), 'R test statistic')
      call near(result%p_value, expected3(2), 'R test p-value')
      if (result%exact) call near(result%point_probability, expected3(3), 'R test point mass')
   end subroutine check_reference3

   subroutine check_reference2()
      read(unit, *) expected2
      call assert(result%status == rank_ok, 'R test status')
      call near(result%statistic, expected2(1), 'R test statistic')
      call near(result%p_value, expected2(2), 'R test p-value')
   end subroutine check_reference2
end program test_exact_rank
