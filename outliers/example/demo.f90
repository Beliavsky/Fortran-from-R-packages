! SPDX-License-Identifier: GPL-2.0-or-later
program demo
   use outliers_api, only : dp, outlier, rm_outlier, scores, grubbs_test, dixon_test, &
      cochran_test, qdixon, qgrubbs, qcochran, outlier_test_result
   implicit none
   real(dp), parameter :: x(6) = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, 20.0_dp]
   type(outlier_test_result) :: test
   print '(a,f8.3)', 'Selected outlier: ', outlier(x)
   print '(a,*(f8.3))', 'Remaining sample: ', rm_outlier(x)
   print '(a,*(f8.3))', 'Normal scores: ', scores(x)
   test = grubbs_test(x)
   if (test%status /= 0) error stop 'Grubbs failed'
   print '(a,2f12.6)', 'Grubbs G and p: ', test%statistic, test%p_value
   test = dixon_test(x)
   if (test%status /= 0) error stop 'Dixon failed'
   print '(a,2f12.6)', 'Dixon Q and p: ', test%statistic, test%p_value
   test = cochran_test([1.0_dp, 2.0_dp, 9.0_dp], [5.0_dp, 5.0_dp, 5.0_dp])
   if (test%status /= 0) error stop 'Cochran failed'
   print '(a,2f12.6)', 'Cochran C and p: ', test%statistic, test%p_value
   print '(a,3f12.6)', 'Critical values: ', qdixon(0.05_dp, 6), qgrubbs(0.95_dp, 6), qcochran(0.95_dp, 5.0_dp, 3)
end program demo
