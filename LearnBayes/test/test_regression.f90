program test_regression
   use learnbayes
   implicit none

   real(dp), parameter :: chirps(15) = [ &
      20.0_dp, 16.0_dp, 19.8_dp, 18.4_dp, 17.1_dp, 15.5_dp, 14.7_dp, 17.1_dp, &
      15.4_dp, 16.2_dp, 15.0_dp, 17.2_dp, 16.0_dp, 17.0_dp, 14.1_dp]
   real(dp), parameter :: temp(15) = [ &
      88.6_dp, 71.6_dp, 93.3_dp, 84.3_dp, 80.6_dp, 75.2_dp, 69.7_dp, 82.0_dp, &
      69.4_dp, 83.3_dp, 78.6_dp, 82.6_dp, 80.6_dp, 83.5_dp, 76.3_dp]
   real(dp) :: x(15, 2)
   real(dp) :: xnew(2, 2)
   real(dp), allocatable :: expected(:, :)
   real(dp), allocatable :: prediction(:, :)
   real(dp) :: beta_mean(2)
   real(dp) :: residual_prob(15)
   type(blinreg_result) :: fit
   type(probit_result) :: pfit
   type(model_selection_result) :: selection
   type(rng_state) :: rng
   integer :: ybin(10)
   real(dp) :: xp(10, 2)
   real(dp) :: prior_beta(2)
   real(dp) :: prior_precision(2, 2)
   integer :: i
   integer :: info

   call rng_seed(rng, 912345_i8)
   x(:, 1) = 1.0_dp
   x(:, 2) = chirps
   call blinreg(rng, temp, x, 6000, fit, info)
   call assert_true(info == 0, 'blinreg status')
   beta_mean = sum(fit%beta, dim=1)/real(size(fit%beta, 1), dp)
   call assert_close(beta_mean(1), 26.01204318_dp, 0.45_dp, 'blinreg intercept mean')
   call assert_close(beta_mean(2), 3.24416574_dp, 0.03_dp, 'blinreg slope mean')
   call assert_true(all(fit%sigma > 0.0_dp), 'blinreg positive sigma')

   xnew = reshape([1.0_dp, 1.0_dp, 16.0_dp, 18.0_dp], [2, 2])
   allocate(expected(size(fit%beta, 1), 2), prediction(size(fit%beta, 1), 2))
   call blinregexpected(xnew, fit, expected)
   call blinregpred(rng, xnew, fit, prediction)
   call assert_true(all(abs(expected) < huge(1.0_dp)), 'blinreg expected finite')
   call assert_true(all(abs(prediction) < huge(1.0_dp)), 'blinreg prediction finite')

   call bayesresiduals(temp, x, fit, 2.0_dp, residual_prob, info)
   call assert_true(info == 0, 'bayesresiduals status')
   call assert_true(all(residual_prob >= 0.0_dp .and. residual_prob <= 1.0_dp), 'bayesresidual probabilities')

   ybin = [0, 1, 0, 0, 0, 1, 1, 1, 1, 1]
   do i = 1, 10
      xp(i, 1) = 1.0_dp
      xp(i, 2) = real(i, dp)
   end do
   prior_beta = 0.0_dp
   prior_precision = 0.0_dp
   prior_precision(1, 1) = 0.5_dp
   prior_precision(2, 2) = 10.0_dp
   call bayes_probit(rng, ybin, xp, 500, pfit, info, prior_beta, prior_precision)
   call assert_true(info == 0, 'bayes.probit status')
   call assert_true(pfit%has_log_marginal, 'bayes.probit marginal flag')
   call assert_true(abs(pfit%log_marginal) < huge(1.0_dp), 'bayes.probit finite marginal')

   call bayes_model_selection(temp, x, 100.0_dp, selection, .true., info)
   call assert_true(info == 0, 'model selection status')
   call assert_close(sum(selection%probability), 1.0_dp, 3.0e-13_dp, 'model probabilities sum')
   call assert_true(size(selection%probability) == 2, 'two candidate models')

   print '(a)', 'test_regression: PASS'

contains

   subroutine assert_close(actual, reference, tolerance, label)
      real(dp), intent(in) :: actual !! Computed scalar value under test.
      real(dp), intent(in) :: reference !! Reference value used for the comparison.
      real(dp), intent(in) :: tolerance !! Maximum permitted absolute difference.
      character(len=*), intent(in) :: label !! Short test label reported on failure.

      if (abs(actual - reference) > tolerance) then
         write (*, '(a,2es24.15)') trim(label)//' failed: ', actual, reference
         error stop 1
      end if
   end subroutine assert_close

   subroutine assert_true(condition, label)
      logical, intent(in) :: condition !! Boolean condition that must be true for the test to pass.
      character(len=*), intent(in) :: label !! Short test label reported when condition is false.

      if (.not. condition) then
         write (*, '(a)') trim(label)//' failed'
         error stop 1
      end if
   end subroutine assert_true

end program test_regression
