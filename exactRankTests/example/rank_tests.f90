! SPDX-License-Identifier: GPL-2.0-or-later
program rank_tests_example
   use exact_rank_tests
   implicit none
   type(rank_test_result) :: test
   type(rank_distribution) :: law
   real(dp), allocatable :: scores(:)
   integer :: status, i
   call wilcox_exact([1.0_dp, 2.0_dp, 2.0_dp, 3.0_dp, 5.0_dp], test, &
      y = [0.0_dp, 2.0_dp, 4.0_dp, 4.0_dp], exact = .true., conf_level = 0.8_dp)
   if (test%status /= rank_ok) error stop 'Wilcoxon calculation failed'
   print '(a,f7.3,a,f9.6)', 'Wilcoxon W = ', test%statistic, ', exact p = ', test%p_value
   print '(a,2f8.3)', '80% location interval: ', test%conf_int
   call cscores([1.0_dp, 2.0_dp, 2.0_dp, 4.0_dp], scores, status, 'Wilcoxon')
   if (status /= rank_ok) error stop 'Scoring failed'
   call permutation_distribution(scores, 2, law, status)
   if (status /= rank_ok) error stop 'Distribution construction failed'
   print '(a)', 'Rank-sum support and exact mass:'
   do i = 1, size(law%support)
      print '(f7.2,f12.7)', law%support(i), law%probability(i)
   end do
end program rank_tests_example
