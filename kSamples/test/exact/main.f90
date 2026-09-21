program test_exact
   use ksamples, only : dp, ad_result, qn_result, jt_result, steel_result, contingency_result
   use ksamples, only : sample_block, contingency_block, combined_result, steel_confint_result
   use ksamples, only : ad_test, ad_test_combined, qn_test, qn_test_combined
   use ksamples, only : jt_test, steel_test, steel_confint, contingency2xt, contingency2xt_comb
   implicit none
   real(dp) :: z(12), jtx(6), steelx(6)
   integer :: ns3(3), ns2(3), avec(2), bvec(2)
   type(ad_result) :: ad
   type(qn_result) :: qn
   type(jt_result) :: jt
   type(steel_result) :: steel
   type(contingency_result) :: table
   type(sample_block) :: blocks(2)
   type(contingency_block) :: tables(2)
   type(combined_result) :: combined
   type(ad_result) :: ad_combined
   type(steel_confint_result) :: interval

   z = [0.824_dp, 0.216_dp, 0.538_dp, 0.685_dp, &
        0.448_dp, 0.348_dp, 0.443_dp, 0.722_dp, &
        0.403_dp, 0.268_dp, 0.440_dp, 0.087_dp]
   ns3 = [4, 4, 4]
   call qn_test(z, ns3, qn, score='KW', method='exact', nsim=100000, return_dist=.false.)
   call assert_close(qn%statistic, 3.576923076923_dp, 1.0e-10_dp, 'QN statistic')
   call assert_close(qn%randomization_p, 0.172987012987_dp, 1.0e-10_dp, 'QN exact p-value')
   call assert_true(trim(qn%method) == 'exact', 'QN exact method')

   call ad_test(z, ns3, method='exact', nsim=100000, keep_null=.false., result=ad)
   call assert_close(ad%statistic(1), 2.636652236652_dp, 2.0e-10_dp, 'AD version 1 statistic')
   call assert_close(ad%statistic(2), 2.624039800251_dp, 2.0e-10_dp, 'AD version 2 statistic')
   call assert_close(ad%randomization_p(1), 0.209235209235_dp, 2.0e-10_dp, 'AD version 1 exact p-value')
   call assert_close(ad%randomization_p(2), 0.217662337662_dp, 2.0e-10_dp, 'AD version 2 exact p-value')

   jtx = [1.0_dp, 2.0_dp, 1.5_dp, 2.1_dp, 1.9_dp, 3.1_dp]
   ns2 = [2, 2, 2]
   call jt_test(jtx, ns2, jt, method='exact', nsim=90, return_dist=.true.)
   call assert_close(jt%statistic, 9.0_dp, 1.0e-14_dp, 'JT statistic')
   call assert_close(jt%mean, 6.0_dp, 1.0e-14_dp, 'JT mean')
   call assert_close(jt%sigma, 2.516611478424_dp, 2.0e-10_dp, 'JT sigma')
   call assert_close(jt%randomization_p, 1.0_dp/6.0_dp, 2.0e-14_dp, 'JT exact p-value')
   call assert_true(size(jt%null_dist) == 90, 'JT exact null size')

   steelx = [1.0_dp, 4.0_dp, 2.0_dp, 5.0_dp, 3.0_dp, 6.0_dp]
   call steel_test(steelx, ns2, alternative='greater', method='exact', nsim=90, &
      keep_null=.true., result=steel)
   call assert_close(steel%statistic, 0.774596669241_dp, 2.0e-10_dp, 'Steel statistic')
   call assert_close(steel%randomization_p, 44.0_dp/90.0_dp, 2.0e-14_dp, 'Steel exact p-value')
   call assert_true(size(steel%null_dist) == 90, 'Steel exact null size')

   avec = [2, 0]
   bvec = [0, 2]
   call contingency2xt(avec, bvec, table, method='exact', nsim=100, return_dist=.true.)
   call assert_close(table%statistic, 3.0_dp, 1.0e-14_dp, 'contingency statistic')
   call assert_close(table%randomization_p, 1.0_dp/3.0_dp, 2.0e-14_dp, 'contingency exact p-value')
   call assert_close(sum(table%probability), 1.0_dp, 2.0e-14_dp, 'contingency probability sum')


   blocks(1)%x = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp]
   blocks(1)%ns = [2, 2]
   blocks(2)%x = [1.0_dp, 4.0_dp, 2.0_dp, 3.0_dp]
   blocks(2)%ns = [2, 2]
   call qn_test_combined(blocks, combined, score='KW', method='exact', nsim=100, return_dist=.true.)
   call assert_true(trim(combined%method) == 'exact', 'combined QN exact method')
   call assert_true(combined%randomization_p >= 0.0_dp .and. combined%randomization_p <= 1.0_dp, &
      'combined QN exact p-value range')

   call ad_test_combined(blocks, method='exact', nsim=100, keep_null=.true., result=ad_combined)
   call assert_true(trim(ad_combined%method) == 'exact', 'combined AD exact method')
   call assert_true(size(ad_combined%null_dist, 1) == 36, 'combined AD exact null size')

   tables(1)%avec = [2, 0]
   tables(1)%bvec = [0, 2]
   tables(2)%avec = [2, 0]
   tables(2)%bvec = [0, 2]
   call contingency2xt_comb(tables, combined, method='exact', nsim=100, return_dist=.true.)
   call assert_close(combined%statistic, 6.0_dp, 1.0e-14_dp, 'combined contingency statistic')
   call assert_close(combined%randomization_p, 1.0_dp/9.0_dp, 2.0e-14_dp, &
      'combined contingency exact p-value')

   call steel_confint(steelx, ns2, conf_level=0.95_dp, alternative='two-sided', &
      method='asymptotic', nsim=100, result=interval)
   call assert_true(trim(interval%method) == 'asymptotic', 'Steel confidence interval method')
   call assert_true(size(interval%lower_conservative) == 2, 'Steel confidence interval size')
   call assert_true(all(interval%lower_conservative <= interval%upper_conservative), &
      'Steel confidence interval ordering')

   print '(a)', 'exact kSamples tests passed'

contains

   subroutine assert_close(actual, expected, tolerance, label)
      real(dp), intent(in) :: actual !! Computed value being checked.
      real(dp), intent(in) :: expected !! Reference value for the check.
      real(dp), intent(in) :: tolerance !! Maximum allowed absolute difference.
      character(len=*), intent(in) :: label !! Human-readable test label shown on failure.
      if (abs(actual - expected) > tolerance) then
         write (*, '(a,2(1x,es24.16))') trim(label)//' failed:', actual, expected
         error stop 1
      end if
   end subroutine assert_close

   subroutine assert_true(condition, label)
      logical, intent(in) :: condition !! Logical condition that must hold.
      character(len=*), intent(in) :: label !! Human-readable test label shown on failure.
      if (.not. condition) then
         write (*, '(a)') trim(label)//' failed'
         error stop 1
      end if
   end subroutine assert_true

end program test_exact
