program test_scale_tests
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_is_nan
   use r_stats, only: dp, fligner_test, chisq_test_result_t
   use stats_test_assertions, only: assert_close, assert_integer, assert_true
   implicit none
   real(dp), parameter :: x(15) = [4.0_dp, 5.0_dp, 6.0_dp, 5.0_dp, 7.0_dp, &
      8.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, 6.0_dp, 6.0_dp, 7.0_dp, 8.0_dp, 9.0_dp, 11.0_dp]
   integer, parameter :: g(15) = [-7, -7, -7, -7, -7, 0, 0, 0, 0, 0, 1000, 1000, 1000, 1000, 1000]
   type(chisq_test_result_t) :: test

   call check(fligner_test(x, g), 1.1631155025200568_dp, 2, 0.55902686314840166_dp, 'ties')
   call check(fligner_test(3.0_dp*x + 100.0_dp, g), &
              1.1631155025200568_dp, 2, 0.55902686314840166_dp, 'location-scale invariance')
   call check(fligner_test(x(15:1:-1), g(15:1:-1)), &
              1.1631155025200568_dp, 2, 0.55902686314840166_dp, 'permutation invariance')
   test = fligner_test([1.0_dp, 2.0_dp, 4.0_dp, 10.0_dp, 12.0_dp, 18.0_dp, 25.0_dp, &
                       -4.0_dp, -3.0_dp, 0.0_dp, 7.0_dp, 9.0_dp], [1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3])
   call check(test, 3.7050292256473902_dp, 2, 0.15684227242955162_dp, 'unequal group sizes')
   test = fligner_test([1.0_dp, 3.0_dp, 6.0_dp, 9.0_dp, 11.0_dp], [1, 2, 2, 2, 2])
   call check(test, 1.6589480870642404_dp, 1, 0.19774541088880235_dp, 'singleton group')
   test = fligner_test([1.0_dp, 1.0_dp, 1.0_dp, 4.0_dp, 4.0_dp, 4.0_dp], [1, 1, 1, 2, 2, 2])
   call assert_true(ieee_is_nan(test%statistic) .and. ieee_is_nan(test%p_value), 'zero score variance')
   print *, 'test_scale_tests: PASS'

contains

   subroutine check(test, statistic, df, probability, label)
      !! Compares a scale test with deterministic R output, explicitly rejecting NaNs.
      type(chisq_test_result_t), intent(in) :: test !! Computed test result.
      real(dp), intent(in) :: statistic !! Reference statistic.
      integer, intent(in) :: df !! Reference degrees of freedom.
      real(dp), intent(in) :: probability !! Reference p-value.
      character(len=*), intent(in) :: label !! Description of the reference case.

      call assert_true(ieee_is_finite(test%statistic) .and. ieee_is_finite(test%p_value), label//' finite')
      call assert_close(test%statistic, statistic, label//' statistic', 2.0e-10_dp)
      call assert_close(test%p_value, probability, label//' probability', 2.0e-10_dp)
      call assert_integer(test%parameter, df, label//' df')
   end subroutine check

end program test_scale_tests
