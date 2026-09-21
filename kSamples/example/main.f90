program rank_tests_example
   use ksamples, only : dp, qn_result, jt_result, qn_test, jt_test
   implicit none
   real(dp) :: z(12), trend(6)
   integer :: ns_qn(3), ns_jt(3)
   type(qn_result) :: qn
   type(jt_result) :: jt

   z = [0.824_dp, 0.216_dp, 0.538_dp, 0.685_dp, &
        0.448_dp, 0.348_dp, 0.443_dp, 0.722_dp, &
        0.403_dp, 0.268_dp, 0.440_dp, 0.087_dp]
   ns_qn = [4, 4, 4]
   call qn_test(z, ns_qn, qn, score='KW', method='exact', nsim=100000)

   trend = [1.0_dp, 2.0_dp, 1.5_dp, 2.1_dp, 1.9_dp, 3.1_dp]
   ns_jt = [2, 2, 2]
   call jt_test(trend, ns_jt, jt, method='exact', nsim=90)

   write (*, '(a,f12.7,a,f12.7)') 'QN statistic = ', qn%statistic, ', exact p = ', qn%randomization_p
   write (*, '(a,f12.7,a,f12.7)') 'JT statistic = ', jt%statistic, ', exact p = ', jt%randomization_p
end program rank_tests_example
