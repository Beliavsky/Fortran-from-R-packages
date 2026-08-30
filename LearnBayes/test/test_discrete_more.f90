program test_discrete_more
   use learnbayes
   implicit none

   type(mixture_beta_result) :: beta_mix
   type(mixture_beta_result) :: gamma_mix
   type(mixture_normal_result) :: normal_mix
   type(discrete_summary) :: summary
   real(dp) :: bf
   real(dp) :: post
   real(dp) :: post_disc(3)
   real(dp) :: pred(4)
   real(dp) :: prior2(2, 3)
   real(dp) :: hist_value(4)

   call pbetat(0.5_dp, 0.25_dp, 1.0_dp, 1.0_dp, 2, 1, bf, post)
   call assert_close(bf, 1.5_dp, 2.0e-14_dp, 'beta point-null Bayes factor')
   call assert_close(post, 1.0_dp/3.0_dp, 2.0e-14_dp, 'beta point-null posterior')

   call pdisc([0.25_dp, 0.5_dp, 0.75_dp], [1.0_dp, 1.0_dp, 1.0_dp]/3.0_dp, 2, 1, post_disc)
   call assert_close(sum(post_disc), 1.0_dp, 2.0e-14_dp, 'discrete posterior normalization')
   call assert_true(post_disc(3) > post_disc(1), 'discrete posterior ordering')

   call pdiscp([0.25_dp, 0.5_dp, 0.75_dp], post_disc, 3, [0, 1, 2, 3], pred)
   call assert_close(sum(pred), 1.0_dp, 2.0e-14_dp, 'discrete predictive normalization')

   call normal_normal_mix([0.4_dp, 0.6_dp], reshape([-1.0_dp, 1.0_dp, 1.0_dp, 4.0_dp], [2, 2]), &
      0.5_dp, 1.0_dp, normal_mix)
   call assert_close(sum(normal_mix%probs), 1.0_dp, 2.0e-14_dp, 'normal mixture weights')
   call assert_true(all(normal_mix%par(:, 2) > 0.0_dp), 'normal mixture posterior variances')

   call binomial_beta_mix([0.3_dp, 0.7_dp], reshape([1.0_dp, 5.0_dp, 1.0_dp, 2.0_dp], [2, 2]), 4, 2, beta_mix)
   call assert_close(sum(beta_mix%probs), 1.0_dp, 2.0e-14_dp, 'beta mixture weights')
   call assert_true(all(beta_mix%par > 0.0_dp), 'beta mixture posterior parameters')

   call poisson_gamma_mix([0.5_dp, 0.5_dp], reshape([2.0_dp, 5.0_dp, 1.0_dp, 2.0_dp], [2, 2]), &
      [2, 1], [1.0_dp, 2.0_dp], gamma_mix)
   call assert_close(sum(gamma_mix%probs), 1.0_dp, 2.0e-14_dp, 'gamma mixture weights')
   call assert_true(all(gamma_mix%par > 0.0_dp), 'gamma mixture posterior parameters')

   call prior_two_parameters([1.0_dp, 2.0_dp], [3.0_dp, 4.0_dp, 5.0_dp], prior2)
   call assert_close(sum(prior2), 1.0_dp, 2.0e-14_dp, 'two-parameter uniform prior')

   call histprior([0.1_dp, 0.3_dp, 0.6_dp, 0.9_dp], [0.25_dp, 0.75_dp], [0.4_dp, 0.6_dp], hist_value)
   call assert_true(all(hist_value >= 0.0_dp), 'histogram prior nonnegative')

   summary = summarize_discrete([1.0_dp, 2.0_dp, 3.0_dp], [0.2_dp, 0.5_dp, 0.3_dp], 0.7_dp)
   call assert_close(summary%mean, 2.1_dp, 2.0e-14_dp, 'discrete summary mean')
   call assert_true(summary%coverage >= 0.7_dp, 'discrete summary HPD coverage')

   print '(a)', 'test_discrete_more: PASS'

contains

   subroutine assert_close(actual, expected, tol, label)
      real(dp), intent(in) :: actual !! Computed scalar value under test.
      real(dp), intent(in) :: expected !! Reference scalar value expected from the translated algorithm.
      real(dp), intent(in) :: tol !! Maximum permitted absolute error.
      character(len=*), intent(in) :: label !! Short test label reported if the comparison fails.

      if (abs(actual - expected) > tol) then
         write (*, '(a,2es24.15)') trim(label)//' failed: ', actual, expected
         error stop 1
      end if
   end subroutine assert_close

   subroutine assert_true(condition, label)
      logical, intent(in) :: condition !! Boolean condition that must be true for the test to pass.
      character(len=*), intent(in) :: label !! Short test label reported if the condition is false.

      if (.not. condition) then
         write (*, '(a)') trim(label)//' failed'
         error stop 1
      end if
   end subroutine assert_true

end program test_discrete_more
