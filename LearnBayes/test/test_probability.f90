program test_probability
   use learnbayes
   implicit none

   real(dp) :: cov2(2, 2)
   real(dp) :: zero2(2)
   real(dp) :: probs(6)
   real(dp) :: ab(2)
   real(dp) :: mu
   real(dp) :: sigma
   real(dp) :: bf
   real(dp) :: prior_odds
   real(dp) :: post_odds
   real(dp) :: post_h
   real(dp) :: ytab(2, 2)
   real(dp) :: atab(2, 2)
   real(dp) :: data_bb(3, 2)
   integer :: s
   integer :: info

   cov2 = 0.0_dp
   cov2(1, 1) = 1.0_dp
   cov2(2, 2) = 1.0_dp
   zero2 = 0.0_dp
   call assert_close(dmnorm(zero2, zero2, cov2), 0.15915494309189535_dp, 2.0e-14_dp, 'dmnorm')
   call assert_close(regularized_beta(0.5_dp, 2.0_dp, 2.0_dp), 0.5_dp, 2.0e-13_dp, 'beta cdf')
   call assert_close(normal_quantile(0.975_dp), 1.959963984540054_dp, 3.0e-9_dp, 'normal quantile')

   do s = 0, 5
      probs(s + 1) = pbetap(1.0_dp, 1.0_dp, 5, s)
   end do
   call assert_true(maxval(abs(probs - 1.0_dp/6.0_dp)) < 2.0e-14_dp, 'uniform beta-binomial predictive')

   call normal_select(0.025_dp, -1.959963984540054_dp, 0.975_dp, 1.959963984540054_dp, mu, sigma)
   call assert_close(mu, 0.0_dp, 1.0e-8_dp, 'normal.select mu')
   call assert_close(sigma, 1.0_dp, 1.0e-8_dp, 'normal.select sigma')

   call beta_select(0.25_dp, 0.25_dp, 0.75_dp, 0.75_dp, ab, info)
   call assert_true(info == 0, 'beta.select status')
   call assert_true(all(ab > 0.0_dp), 'beta.select positive shapes')

   ytab = reshape([3.0_dp, 2.0_dp, 1.0_dp, 4.0_dp], [2, 2])
   atab = 1.0_dp
   call assert_close(ctable(ytab, atab), 1.7769230769230846_dp, 5.0e-13_dp, 'ctable')

   data_bb = reshape([2.0_dp, 4.0_dp, 1.0_dp, 5.0_dp, 6.0_dp, 4.0_dp], [3, 2])
   call assert_close(betabinexch([0.2_dp, log(5.0_dp)], data_bb), -13.10191433834591_dp, 3.0e-13_dp, 'betabinexch')

   call mnormt_onesided(0.0_dp, 0.0_dp, 2.0_dp, 0.5_dp, 10, 1.0_dp, bf, prior_odds, post_odds, post_h)
   call assert_close(bf, 0.06289681978970314_dp, 2.0e-12_dp, 'mnormt.onesided BF')
   call assert_close(post_h, 0.059174906367814174_dp, 2.0e-12_dp, 'mnormt.onesided probability')

   call assert_close(logctablepost([0.4_dp, -0.2_dp], [3.0_dp, 7.0_dp, 4.0_dp, 6.0_dp]), &
      -13.887519045420982_dp, 2.0e-13_dp, 'logctablepost')

   print '(a)', 'test_probability: PASS'

contains

   subroutine assert_close(actual, expected, tol, label)
      real(dp), intent(in) :: actual !! Computed scalar value under test.
      real(dp), intent(in) :: expected !! Reference scalar value expected from the upstream formula.
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

end program test_probability
