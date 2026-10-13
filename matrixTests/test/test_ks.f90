! SPDX-License-Identifier: GPL-2.0-only
program test_ks
   !! R fixtures were evaluated with matrixTests 0.2.3.1 under R 4.5.0.
   !! Independent exhaustive tests enumerate every allocation of pooled labels.
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, ieee_positive_inf, ieee_is_nan
   use r_kinds, only : dp
   use matrixtests_ks, only : ks_test, ks_success, ks_invalid_sample, ks_invalid_alternative, ks_resource_limit
   implicit none
   character(9), parameter :: alternatives(3) = [character(9) :: "two.sided", "less", "greater"]
   integer :: failures, n, i, status
   real(dp) :: d, p, nan, infinity, reference, empty(0)
   real(dp) :: lower(500), upper(500), large(5001)
   logical :: used_exact

   failures = 0
   nan = ieee_value(0.0_dp, ieee_quiet_nan)
   infinity = ieee_value(0.0_dp, ieee_positive_inf)
   call r_fixture([1.0_dp, 4.0_dp, 6.0_dp], [2.0_dp, 3.0_dp, 5.0_dp, 7.0_dp], &
      [1.0_dp / 3.0_dp, 1.0_dp / 6.0_dp, 1.0_dp / 3.0_dp], &
      [0.97142857142857131_dp, 0.80000000000000004_dp, 0.59999999999999998_dp], &
      [0.99116359639132134_dp, 0.90915644287671304_dp, 0.68321042267494825_dp])
   call r_fixture([1.0_dp, 1.0_dp, 2.0_dp, 4.0_dp], [1.0_dp, 2.0_dp, 2.0_dp, 3.0_dp, 4.0_dp], &
      [0.3_dp, 0.05_dp, 0.3_dp], &
      [0.85714285714285721_dp, 0.82539682539682535_dp, 0.50793650793650791_dp], &
      [0.98826107764352444_dp, 0.98895038929392232_dp, 0.67032004603563933_dp])
   call r_fixture([0.0_dp, 0.0_dp, 0.0_dp], [1.0_dp, 1.0_dp, 2.0_dp, 2.0_dp], &
      [1.0_dp, 0.0_dp, 1.0_dp], &
      [0.028571428571428571_dp, 1.0_dp, 0.028571428571428571_dp], &
      [0.064864268738852426_dp, 1.0_dp, 0.03243324089479549_dp])
   call r_fixture([1.0_dp, 2.0_dp, 2.0_dp, 4.0_dp, 8.0_dp], &
      [1.0_dp, 1.0_dp, 2.0_dp, 3.0_dp, 3.0_dp, 8.0_dp, 9.0_dp], &
      [0.17142857142857149_dp, 0.11428571428571421_dp, 0.17142857142857149_dp], &
      [1.0_dp, 0.72979797979797978_dp, 0.63636363636363635_dp], &
      [0.99999519506318424_dp, 0.92663968695101029_dp, 0.84246044161677125_dp])

   ! Every allocation is tested, with both unequal sample sizes and tied values.
   do n = 1, 8
      call exhaustive_labels([real(dp) :: 5, 1, 9, 3, 6, 2, 7, 4, 8], n)
      call exhaustive_labels([real(dp) :: 1, 0, 2, 0, 1, 0, 2, 1, 2], n)
   end do
   call exhaustive_labels([-infinity, 0.0_dp, -infinity, 2.0_dp, -infinity, 1.0_dp, infinity], 4)
   call exhaustive_labels([real(dp) :: 1, 1, 1, 1, 1, 1, 1], 3)
   call ks_test([-infinity, -infinity, 0.0_dp, 2.0_dp], [-infinity, 1.0_dp, infinity], &
      d, p, status, exact = .true.)
   call require(status == ks_success, "infinities accepted")
   call close(d, 5.0_dp / 12.0_dp, 1.0e-15_dp, "infinite ties statistic")
   ! R 4.5.0 gives 31/35 here because diff(-Inf,-Inf) is NaN. Exhaustive
   ! conditional-label enumeration gives 29/35, matching finite rank equivalents.
   call close(p, 29.0_dp / 35.0_dp, 2.0e-15_dp, "infinite ties exact-label probability")

   ! A directly evaluated exact survival tail remains nonzero below 1e-298.
   lower = [(real(i, dp), i = 1, 500)]
   upper = lower + 500.0_dp
   ! 70-decimal-digit independent Jacobi-series references. R 4.5.0's
   ! tol=1e-6 small-z truncation returns 0.39273380987443274 and
   ! 0.28095995364725834, respectively, rather than these converged values.
   call ks_test(lower(:50), lower(:50) + 9.0_dp, d, p, status, exact = .false.)
   call close(p, 0.3927307079406543739271028715_dp, 5.0e-15_dp, "converged limit at scaled statistic .9")
   call ks_test(lower(:100), lower(:100) + 14.0_dp, d, p, status, exact = .false.)
   call close(p, 0.2809295474198877036889821833_dp, 5.0e-15_dp, "converged limit near scaled statistic 1")
   reference = 2.0_dp
   do i = 1, 500
      reference = reference * real(i, dp) / real(500 + i, dp)
   end do
   call ks_test(lower, upper, d, p, status, exact = .true., used_exact = used_exact)
   call require(status == ks_success .and. used_exact, "extreme exact tail succeeds")
   call require(p > 0.0_dp, "extreme exact tail avoids cancellation")
   call close(p / reference, 1.0_dp, 2.0e-13_dp, "extreme exact tail ratio")
   call ks_test(lower, upper, d, p, status, exact = .false.)
   call close(p / (2.0_dp * exp(-500.0_dp)), 1.0_dp, 2.0e-13_dp, "extreme two-sided asymptotic tail")
   call ks_test(lower, upper, d, p, status, alternative = "greater", exact = .false.)
   call close(p / exp(-500.0_dp), 1.0_dp, 2.0e-13_dp, "extreme one-sided asymptotic tail")

   call ks_test(lower(:99), upper(:101), d, p, status, used_exact = used_exact)
   call require(status == ks_success .and. used_exact, "automatic exact below product 10000")
   call ks_test(lower(:100), upper(:100), d, p, status, used_exact = used_exact)
   call require(status == ks_success .and. .not. used_exact, "automatic asymptotic at product 10000")
   call ks_test(lower(:100), upper(:100), d, p, status, exact = .true., used_exact = used_exact)
   call require(status == ks_success .and. used_exact, "explicit exact overrides automatic selection")
   call ks_test([1.0_dp], [2.0_dp], d, p, status, alternative = "g", exact = .true.)
   call close(p, 0.5_dp, 1.0e-15_dp, "R-style alternative prefix")
   call ks_test(empty, [1.0_dp], d, p, status)
   call require(status == ks_invalid_sample .and. ieee_is_nan(p), "empty sample rejected")
   call ks_test([nan], [1.0_dp], d, p, status)
   call require(status == ks_invalid_sample .and. ieee_is_nan(p), "kernel NaN rejected")
   call ks_test([1.0_dp], [2.0_dp], d, p, status, alternative = "invalid")
   call require(status == ks_invalid_alternative .and. ieee_is_nan(p), "invalid alternative rejected")
   call ks_test([1.0_dp], [2.0_dp], d, p, status, alternative = "")
   call require(status == ks_invalid_alternative, "empty alternative rejected")
   large = 1.0_dp
   call ks_test(large, large, d, p, status, exact = .true., used_exact = used_exact)
   call require(status == ks_resource_limit .and. ieee_is_nan(p), "oversized exact request reported")
   call require(.not. used_exact, "failed exact request never claims an exact result")

   if (failures /= 0) error stop "Kolmogorov-Smirnov tests failed"
   print '(a)', "Kolmogorov-Smirnov R fixtures and exhaustive-label tests passed."

contains

   subroutine require(condition, label)
      logical, intent(in) :: condition !! Assertion that must hold for a successful test.
      character(*), intent(in) :: label !! Human-readable identification of the tested property.

      if (condition) return
      failures = failures + 1
      if (failures <= 12) print '(a)', "FAIL: " // label
   end subroutine require

   subroutine close(actual, expected, tolerance, label)
      real(dp), intent(in) :: actual !! Computed numerical value to validate.
      real(dp), intent(in) :: expected !! Independent reference value.
      real(dp), intent(in) :: tolerance !! Nonnegative maximum absolute difference.
      character(*), intent(in) :: label !! Human-readable identification of the comparison.

      call require(.not. ieee_is_nan(actual), label // " is not NaN")
      call require(abs(actual - expected) <= tolerance, label)
   end subroutine close

   subroutine r_fixture(x, y, statistics, exact_probabilities, asymptotic_probabilities)
      real(dp), intent(in) :: x(:) !! First sample supplied to the installed R matrixTests oracle.
      real(dp), intent(in) :: y(:) !! Second sample supplied to the installed R matrixTests oracle.
      real(dp), intent(in) :: statistics(3) !! R statistics in two.sided, less, greater order.
      real(dp), intent(in) :: exact_probabilities(3) !! R exact p-values in two.sided, less, greater order.
      real(dp), intent(in) :: asymptotic_probabilities(3) !! R limiting p-values in two.sided, less, greater order.
      real(dp) :: statistic, probability
      integer :: alternative, result_status
      logical :: selected_exact

      do alternative = 1, 3
         call ks_test(x, y, statistic, probability, result_status, alternatives(alternative), .true., selected_exact)
         call require(result_status == ks_success .and. selected_exact, "R fixture exact status")
         call close(statistic, statistics(alternative), 5.0e-15_dp, "R fixture statistic")
         call close(probability, exact_probabilities(alternative), 5.0e-14_dp, "R fixture exact p-value")
         call ks_test(x, y, statistic, probability, result_status, alternatives(alternative), .false., selected_exact)
         call require(result_status == ks_success .and. .not. selected_exact, "R fixture asymptotic status")
         call close(probability, asymptotic_probabilities(alternative), 5.0e-13_dp, "R fixture asymptotic p-value")
      end do
   end subroutine r_fixture

   subroutine exhaustive_labels(pool, nx)
      !! Build an independent integer statistic histogram for all choose(N,nx)
      !! labelled partitions, then compare every partition with the exact kernel.
      !! The reference ECDF is evaluated by COUNT at each observed value, without
      !! sorting, a lattice recursion, or a shared statistic helper.
      real(dp), intent(in) :: pool(:) !! Unsorted pooled sample, at most nine values, with ties or infinities permitted.
      integer, intent(in) :: nx !! Positive number assigned to x, strictly below the pooled sample size.
      integer, allocatable :: histogram(:, :)
      real(dp), allocatable :: x(:), y(:)
      logical, allocatable :: assignment(:)
      integer :: total, ny, mask, k, a, assignments, deviation, thresholds(3), result_status
      real(dp) :: statistic, probability, expected

      total = size(pool)
      ny = total - nx
      allocate(histogram(0:nx * ny, 3), x(nx), y(ny), assignment(total))
      histogram = 0
      assignments = 0
      do mask = 0, 2**total - 1
         if (popcnt(mask) /= nx) cycle
         assignment = [(btest(mask, k - 1), k = 1, total)]
         x = pack(pool, assignment)
         y = pack(pool, .not. assignment)
         thresholds = 0
         do k = 1, total
            deviation = count(x <= pool(k)) * ny - count(y <= pool(k)) * nx
            thresholds(1) = max(thresholds(1), abs(deviation))
            thresholds(2) = max(thresholds(2), -deviation)
            thresholds(3) = max(thresholds(3), deviation)
         end do
         do a = 1, 3
            histogram(thresholds(a), a) = histogram(thresholds(a), a) + 1
         end do
         assignments = assignments + 1
      end do
      do mask = 0, 2**total - 1
         if (popcnt(mask) /= nx) cycle
         assignment = [(btest(mask, k - 1), k = 1, total)]
         x = pack(pool, assignment)
         y = pack(pool, .not. assignment)
         thresholds = 0
         do k = 1, total
            deviation = count(x <= pool(k)) * ny - count(y <= pool(k)) * nx
            thresholds(1) = max(thresholds(1), abs(deviation))
            thresholds(2) = max(thresholds(2), -deviation)
            thresholds(3) = max(thresholds(3), deviation)
         end do
         do a = 1, 3
            expected = real(sum(histogram(thresholds(a):, a)), dp) / real(assignments, dp)
            call ks_test(x, y, statistic, probability, result_status, alternatives(a), .true.)
            call require(result_status == ks_success, "exhaustive exact status")
            call close(statistic, real(thresholds(a), dp) / real(nx * ny, dp), 1.0e-15_dp, &
               "independent ECDF statistic")
            call close(probability, expected, 3.0e-15_dp, "enumerated conditional-label p-value")
         end do
      end do
   end subroutine exhaustive_labels

end program test_ks
