program test_mcnemar
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_is_nan
   use r_stats, only: chisq_test_result_t, dp, mcnemar_test
   use stats_test_assertions, only: assert_close, assert_integer, assert_true
   implicit none
   integer, parameter :: binary(2, 2) = reshape([20, 12, 5, 30], [2, 2])
   integer, parameter :: three(3, 3) = reshape([10, 4, 2, 7, 12, 6, 3, 9, 8], [3, 3])
   real(dp), parameter :: fractional(2, 2) = reshape([20.0_dp, 1.25_dp, 1.5_dp, 30.0_dp], [2, 2])
   type(chisq_test_result_t) :: test

   call check(mcnemar_test(binary), 2.1176470588235294_dp, 1, 0.14561009539686795_dp, 'corrected')
   call check(mcnemar_test(binary, .false.), 2.8823529411764706_dp, 1, 0.089555074413643118_dp, 'uncorrected')
   call check(mcnemar_test(transpose(binary)), 2.1176470588235294_dp, 1, 0.14561009539686795_dp, 'transpose')
   call check(mcnemar_test(reshape([20, 7, 7, 30], [2, 2])), 0.0_dp, 1, 1.0_dp, 'equal counts')
   call check(mcnemar_test(fractional), 0.20454545454545456_dp, 1, 0.65107663407783434_dp, 'fractional')
   call check(mcnemar_test(fractional, .false.), 0.022727272727272728_dp, 1, &
              0.88016845490672546_dp, 'fractional uncorrected')
   call check(mcnemar_test(three), 1.6181818181818182_dp, 3, 0.65527427736577593_dp, 'three categories')
   call check(mcnemar_test(three, .false.), 1.6181818181818182_dp, 3, &
              0.65527427736577593_dp, 'three categories ignores correction')
   test = mcnemar_test(reshape([10, 0, 0, 20], [2, 2]))
   call assert_true(ieee_is_nan(test%statistic), 'zero discordance undefined statistic')
   call assert_true(ieee_is_nan(test%p_value), 'zero discordance undefined probability')
   print *, 'test_mcnemar: PASS'

contains

   subroutine check(test, statistic, df, probability, label)
      !! Verifies finite test results against independently generated R fixtures.
      type(chisq_test_result_t), intent(in) :: test !! Computed McNemar result.
      real(dp), intent(in) :: statistic !! Reference statistic.
      integer, intent(in) :: df !! Reference degrees of freedom.
      real(dp), intent(in) :: probability !! Reference p-value.
      character(len=*), intent(in) :: label !! Name of the reference case.

      call assert_true(ieee_is_finite(test%statistic) .and. ieee_is_finite(test%p_value), label//' finite')
      call assert_close(test%statistic, statistic, label//' statistic', 2.0e-12_dp)
      call assert_integer(test%parameter, df, label//' df')
      call assert_close(test%p_value, probability, label//' probability', 2.0e-12_dp)
   end subroutine check

end program test_mcnemar
