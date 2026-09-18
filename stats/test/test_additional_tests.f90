program test_additional_tests
   use r_stats, only: bartlett_test, binom_test, chisq_test_result_t, dp, &
                      exact_count_test_result_t, friedman_test, poisson_test, &
                      var_test, variance_test_result_t
   use stats_test_assertions, only: assert_close, assert_integer
   implicit none

   real(dp), parameter :: first_sample(6) = [4.2_dp, 5.1_dp, 3.8_dp, 6.0_dp, 5.5_dp, 4.7_dp]
   real(dp), parameter :: second_sample(7) = [2.1_dp, 2.8_dp, 3.5_dp, 2.9_dp, 3.2_dp, 2.6_dp, 3.0_dp]
   real(dp), parameter :: grouped_values(15) = [4.0_dp, 5.0_dp, 6.0_dp, 5.0_dp, 7.0_dp, &
                                                8.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, 6.0_dp, &
                                                6.0_dp, 7.0_dp, 8.0_dp, 9.0_dp, 11.0_dp]
   integer, parameter :: groups(15) = [1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3]
   real(dp), parameter :: blocks(5, 4) = reshape([8.0_dp, 5.0_dp, 9.0_dp, 6.0_dp, 7.0_dp, &
                                                  7.0_dp, 5.0_dp, 8.0_dp, 7.0_dp, 6.0_dp, &
                                                  9.0_dp, 7.0_dp, 10.0_dp, 8.0_dp, 8.0_dp, &
                                                  6.0_dp, 4.0_dp, 7.0_dp, 5.0_dp, 6.0_dp], [5, 4])
   type(chisq_test_result_t) :: chi_square
   type(exact_count_test_result_t) :: exact
   type(variance_test_result_t) :: variance

   exact = binom_test(7, 20, probability=0.2_dp, confidence_level=0.9_dp)
   call assert_integer(exact%status, 0, "binomial status")
   call assert_close(exact%p_value, 0.098221728613468728_dp, "R binomial p-value", 2.0e-12_dp)
   call assert_close(exact%estimate, 0.35_dp, "R binomial estimate", 1.0e-15_dp)
   call assert_close(exact%conf_low, 0.17731091757444906_dp, "R binomial lower interval", 2.0e-10_dp)
   call assert_close(exact%conf_high, 0.55803451131548865_dp, "R binomial upper interval", 2.0e-10_dp)

   exact = poisson_test(12, 4.0_dp, rate=2.0_dp, confidence_level=0.9_dp)
   call assert_integer(exact%status, 0, "Poisson status")
   call assert_close(exact%p_value, 0.15430411301020253_dp, "R Poisson p-value", 2.0e-12_dp)
   call assert_close(exact%estimate, 3.0_dp, "R Poisson estimate", 1.0e-15_dp)
   call assert_close(exact%conf_low, 1.7310531283962769_dp, "R Poisson lower interval", 2.0e-10_dp)
   call assert_close(exact%conf_high, 4.8606423324787542_dp, "R Poisson upper interval", 2.0e-10_dp)

   variance = var_test(first_sample, second_sample, ratio=1.25_dp, confidence_level=0.9_dp)
   call assert_integer(variance%status, 0, "variance test status")
   call assert_close(variance%statistic, 2.6914832535885167_dp, "R variance statistic", 2.0e-12_dp)
   call assert_close(variance%p_value, 0.25997191462696345_dp, "R variance p-value", 2.0e-11_dp)
   call assert_close(variance%conf_low, 0.76682633467711914_dp, "R variance lower interval", 2.0e-9_dp)
   call assert_close(variance%conf_high, 16.654521796662248_dp, "R variance upper interval", 2.0e-8_dp)

   chi_square = bartlett_test(grouped_values, groups)
   call assert_close(chi_square%statistic, 1.13436830595397_dp, "R Bartlett statistic", 2.0e-12_dp)
   call assert_integer(chi_square%parameter, 2, "R Bartlett degrees of freedom")
   call assert_close(chi_square%p_value, 0.56712011595245393_dp, "R Bartlett p-value", 2.0e-12_dp)

   chi_square = friedman_test(blocks)
   call assert_close(chi_square%statistic, 13.5625_dp, "R Friedman statistic", 2.0e-12_dp)
   call assert_integer(chi_square%parameter, 3, "R Friedman degrees of freedom")
   call assert_close(chi_square%p_value, 0.0035654012691661513_dp, "R Friedman p-value", 2.0e-12_dp)

   print *, "test_additional_tests: PASS"

end program test_additional_tests
