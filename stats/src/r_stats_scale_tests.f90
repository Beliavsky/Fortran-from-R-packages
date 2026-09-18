! SPDX-License-Identifier: MIT
module r_stats_scale_tests
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_quiet_nan, ieee_value
   use r_kinds, only: dp
   use r_distributions, only: r_pchisq, r_qnorm
   use r_sorting, only: r_average_ranks
   use r_stats_descriptive, only: median
   use r_stats_types, only: chisq_test_result_t
   implicit none
   private
   public :: fligner_test

contains

   pure function fligner_test(x, group) result(test)
      !! Tests equal group scales using normal scores of ranked absolute deviations from group medians.
      real(dp), intent(in) :: x(:) !! Finite observations, one per group label.
      integer, intent(in) :: group(:) !! Arbitrary integer labels defining the groups.
      type(chisq_test_result_t) :: test !! Chi-square statistic, group count minus one, and p-value.
      real(dp), allocatable :: deviations(:), ranks(:), scores(:), sample(:)
      integer, allocatable :: labels(:), membership(:)
      real(dp) :: center, grand_mean, score_variance, group_mean
      integer :: i, j, n, number_groups, group_size

      n = size(x)
      if (size(group) /= n .or. n < 2) error stop 'fligner_test: incompatible or empty samples'
      if (.not. all(ieee_is_finite(x))) error stop 'fligner_test: observations must be finite'
      allocate (labels(n), membership(n), deviations(n))
      number_groups = 0
      do i = 1, n
         j = 1
         do while (j <= number_groups)
            if (labels(j) == group(i)) exit
            j = j + 1
         end do
         if (j > number_groups) then
            number_groups = j
            labels(j) = group(i)
         end if
         membership(i) = j
      end do
      if (number_groups < 2) error stop 'fligner_test: at least two groups are required'
      do j = 1, number_groups
         sample = pack(x, membership == j)
         center = median(sample)
         where (membership == j) deviations = abs(x - center)
      end do
      call r_average_ranks(deviations, ranks)
      scores = r_qnorm(0.5_dp + ranks/(2.0_dp*real(n + 1, dp)))
      grand_mean = sum(scores)/real(n, dp)
      score_variance = sum((scores - grand_mean)**2)/real(n - 1, dp)
      test%parameter = number_groups - 1
      test%method = 6
      if (maxval(deviations) == minval(deviations) .or. score_variance <= 0.0_dp) then
         test%statistic = ieee_value(0.0_dp, ieee_quiet_nan)
         test%p_value = test%statistic
         return
      end if
      test%statistic = 0.0_dp
      do j = 1, number_groups
         group_size = count(membership == j)
         group_mean = sum(scores, mask=membership == j)/real(group_size, dp)
         test%statistic = test%statistic + real(group_size, dp)*(group_mean - grand_mean)**2
      end do
      test%statistic = test%statistic/score_variance
      test%p_value = r_pchisq(test%statistic, real(test%parameter, dp), lower_tail=.false.)
   end function fligner_test

end module r_stats_scale_tests
