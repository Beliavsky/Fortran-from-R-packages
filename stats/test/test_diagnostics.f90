program test_diagnostics
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_is_nan
   use r_stats, only: box_test, chisq_test_result_t, dp
   use stats_test_assertions, only: assert_close, assert_integer, assert_true
   implicit none
   real(dp), parameter :: x(20) = [2.0_dp, 3.0_dp, 5.0_dp, 4.0_dp, 6.0_dp, &
      8.0_dp, 7.0_dp, 9.0_dp, 10.0_dp, 8.0_dp, 7.0_dp, 6.0_dp, 4.0_dp, 5.0_dp, &
      3.0_dp, 2.0_dp, 4.0_dp, 6.0_dp, 5.0_dp, 7.0_dp]
   type(chisq_test_result_t) :: test

   call check(box_test(x), 9.0170190952998848_dp, 1, 0.0026747725248282839_dp, 'default')
   call check(box_test(x, 5), 17.690666607974357_dp, 5, 0.0033601460888317281_dp, 'BP lag 5')
   call check(box_test(x, 5, fitdf=2), 17.690666607974357_dp, 3, &
              0.00050941889059241685_dp, 'BP fitdf')
   call check(box_test(x, 19), 39.500081920972299_dp, 19, 0.0038039410486058189_dp, 'BP max lag')
   call check(box_test(x, 1, 'Ljung-Box'), 10.440758952452498_dp, 1, &
              0.0012326460997807986_dp, 'LB lag 1')
   call check(box_test(x, 5, 'Ljung-Box'), 22.048009725143586_dp, 5, &
              0.00051270855800866766_dp, 'LB lag 5')
   call check(box_test(x, 5, 'Ljung-Box', 2), 22.048009725143586_dp, 3, &
              6.3747767368393937e-5_dp, 'LB fitdf')
   call check(box_test(x, 19, 'Ljung-Box'), 66.185605943733151_dp, 19, &
              3.9131216711396632e-7_dp, 'LB max lag')
   call check(box_test(2.0_dp*x + 50.0_dp, 5, 'Ljung-Box'), 22.048009725143586_dp, 5, &
              0.00051270855800866766_dp, 'location-scale invariance')
   test = box_test([1.0_dp, 1.0_dp, 1.0_dp])
   call assert_true(ieee_is_nan(test%statistic) .and. ieee_is_nan(test%p_value), 'constant series')
   print *, 'test_diagnostics: PASS'

contains

   subroutine check(test, statistic, df, probability, label)
      !! Verifies finite portmanteau test outputs against R reference values.
      type(chisq_test_result_t), intent(in) :: test !! Computed test result.
      real(dp), intent(in) :: statistic !! Reference statistic.
      integer, intent(in) :: df !! Reference degrees of freedom.
      real(dp), intent(in) :: probability !! Reference probability.
      character(len=*), intent(in) :: label !! Name of the reference case.

      call assert_true(ieee_is_finite(test%statistic) .and. ieee_is_finite(test%p_value), label//' finite')
      call assert_close(test%statistic, statistic, label//' statistic', 2.0e-10_dp)
      call assert_integer(test%parameter, df, label//' df')
      call assert_close(test%p_value, probability, label//' probability', 2.0e-12_dp)
   end subroutine check

end program test_diagnostics
