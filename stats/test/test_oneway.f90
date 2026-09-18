program test_oneway
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_is_nan
   use r_stats, only: dp, oneway_test, oneway_test_result_t
   use stats_test_assertions, only: assert_close, assert_integer, assert_true
   implicit none
   real(dp), parameter :: x(12) = [1.0_dp, 2.0_dp, 4.0_dp, 10.0_dp, 12.0_dp, 18.0_dp, &
                                   25.0_dp, -4.0_dp, -3.0_dp, 0.0_dp, 7.0_dp, 9.0_dp]
   integer, parameter :: g(12) = [-7, -7, -7, 0, 0, 0, 0, 1000, 1000, 1000, 1000, 1000]
   type(oneway_test_result_t) :: test

   call check(oneway_test(x, g), 7.1684092622428315_dp, 2.0_dp, &
      5.2048027158130648_dp, 0.031972298244386403_dp, 'Welch three groups')
   call check(oneway_test(x, g, .true.), 8.704514363885087_dp, 2.0_dp, &
      9.0_dp, 0.0078741804285178806_dp, 'pooled three groups')
   call check(oneway_test(x(:7), g(:7)), 15.909298345693097_dp, 1.0_dp, &
      3.399726075590304_dp, 0.022303825181146802_dp, 'Welch two groups')
   call check(oneway_test(x(:7), g(:7), .true.), 11.73878272581867_dp, 1.0_dp, &
      5.0_dp, 0.018712451578049419_dp, 'pooled two groups')
   call check(oneway_test(2.0_dp*x + 30.0_dp, g), 7.1684092622428315_dp, 2.0_dp, &
      5.2048027158130648_dp, 0.031972298244386403_dp, 'location-scale invariance')
   test = oneway_test([1.0_dp, 1.0_dp, 2.0_dp, 3.0_dp], [1, 1, 2, 2])
   call assert_integer(test%status, 2, 'zero within-group variance status')
   call assert_true(ieee_is_nan(test%p_value), 'undefined Welch probability')
   test = oneway_test(x, g(:11))
   call assert_integer(test%status, 1, 'shape mismatch')
   test = oneway_test([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp], [1, 2, 2, 2])
   call assert_integer(test%status, 1, 'singleton group rejected')
   print *, 'test_oneway: PASS'

contains

   subroutine check(test, statistic, df1, df2, probability, label)
      !! Checks both degrees of freedom and the test statistic against R reference values.
      type(oneway_test_result_t), intent(in) :: test !! Computed ANOVA result.
      real(dp), intent(in) :: statistic !! Reference F statistic.
      real(dp), intent(in) :: df1 !! Reference numerator degrees of freedom.
      real(dp), intent(in) :: df2 !! Reference denominator degrees of freedom.
      real(dp), intent(in) :: probability !! Reference p-value.
      character(len=*), intent(in) :: label !! Reference-case name.

      call assert_integer(test%status, 0, label//' status')
      call assert_true(ieee_is_finite(test%statistic) .and. ieee_is_finite(test%p_value), label//' finite')
      call assert_close(test%statistic, statistic, label//' F', 2.0e-10_dp)
      call assert_close(test%df_numerator, df1, label//' numerator df')
      call assert_close(test%df_denominator, df2, label//' denominator df', 2.0e-10_dp)
      call assert_close(test%p_value, probability, label//' p', 2.0e-10_dp)
   end subroutine check

end program test_oneway
