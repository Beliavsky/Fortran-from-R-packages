! SPDX-License-Identifier: GPL-3.0-or-later
! Deterministic reference values generated using upstream loo 2.10.1 in R.
! CRPS fixtures use set.seed(319), permutations=2, x=matrix((1:24)/4-2,6,4),
! x2 columns seq(4,1,length.out=6), seq(-2,2,length.out=6), seq(0,1,length.out=6),
! seq(2,-1,length.out=6), y=c(.3,1.2,-.5,2.1), log_lik=matrix(-log(1:6),6,4).
! LOO weights are obtained from weights(psis(-log_lik),log=FALSE), with analogous
! joint weights from psis(-log_lik-log_lik[permutation,]). No RNG is used here.
program test_metrics
   use r_kinds, only: dp
   use loo_metrics
   use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
   implicit none
   real(dp) :: x(6, 4), x2(6, 4), y(4), weights(6, 4), joint(6, 4, 2), xc(6, 4)
   real(dp) :: invalid(6, 4)
   integer :: permutations(6, 2), i, j, p, bad_perm(6, 2)
   type(predictive_metric_result) :: out

   x = reshape([(real(i, dp) / 4.0_dp - 2.0_dp, i = 1, 24)], [6, 4])
   do i = 1, 6
      x2(i, 1) = 4.0_dp - 3.0_dp * real(i - 1, dp) / 5.0_dp
      x2(i, 2) = -2.0_dp + 4.0_dp * real(i - 1, dp) / 5.0_dp
      x2(i, 3) = real(i - 1, dp) / 5.0_dp
      x2(i, 4) = 2.0_dp - 3.0_dp * real(i - 1, dp) / 5.0_dp
   end do
   y = [0.3_dp, 1.2_dp, -0.5_dp, 2.1_dp]
   permutations = reshape([2, 3, 4, 6, 1, 5, 4, 3, 6, 1, 2, 5], [6, 2])
   do j = 1, 4
      weights(:, j) = [(real(i, dp) / 21.0_dp, i = 1, 6)]
      do p = 1, 2
         joint(:, j, p) = [(real(i * permutations(i, p), dp), i = 1, 6)]
         joint(:, j, p) = joint(:, j, p) / sum(joint(:, j, p))
      end do
   end do

   out = crps(x, x2, y, permutations)
   call check(out%status == 0, 'crps status')
   call near(out%estimate, -0.351041666666667_dp, 'crps estimate')
   call near(out%se, 0.465690245242359_dp, 'crps standard error')
   call near_vector(out%pointwise, [0.3875_dp, -0.266666666666667_dp, -1.6875_dp, 0.1625_dp], 'crps pointwise')
   out = scrps(x, x2, y, permutations)
   call near(out%estimate, -1.17225369549727_dp, 'scrps estimate')
   call near(out%se, 0.243554182759252_dp, 'scrps standard error')
   out = loo_crps(x, x2, y, weights, joint, permutations)
   call near(out%estimate, -0.376316367791843_dp, 'loo_crps estimate')
   call near(out%se, 0.515248009151816_dp, 'loo_crps standard error')
   call near_vector(out%pointwise, &
      [0.172371659634317_dp, -0.170938818565401_dp, -1.88432665260197_dp, 0.377628340365682_dp], 'loo pointwise')
   out = loo_scrps(x, x2, y, weights, joint, permutations)
   call near(out%estimate, -1.16354527208881_dp, 'loo_scrps estimate')
   call near(out%se, 0.297695580458789_dp, 'loo_scrps standard error')

   out = loo_predictive_metric(x, y, weights)
   call near(out%estimate, 1.475_dp, 'mean absolute error')
   call near(out%se, 0.411495690990445_dp, 'mae standard error')
   out = loo_predictive_metric(x, y, weights, 'mse')
   call near(out%estimate, 2.68361111111111_dp, 'mean squared error')
   call near(out%se, 1.38163751477087_dp, 'mse standard error')
   out = loo_predictive_metric(x, y, weights, 'rmse')
   call near(out%estimate, 1.6381731016932_dp, 'root mean squared error')
   call near(out%se, 0.421700708350914_dp, 'rmse standard error')
   xc = reshape([0.1_dp, 0.2_dp, 0.3_dp, 0.7_dp, 0.8_dp, 0.9_dp, &
      0.9_dp, 0.8_dp, 0.7_dp, 0.2_dp, 0.1_dp, 0.0_dp, &
      0.0_dp, 0.1_dp, 0.2_dp, 0.3_dp, 0.4_dp, 0.5_dp, &
      0.9_dp, 0.8_dp, 0.7_dp, 0.6_dp, 0.5_dp, 0.4_dp], [6, 4])
   out = loo_predictive_metric(xc, [0.0_dp, 0.0_dp, 1.0_dp, 1.0_dp], weights, 'acc')
   call near(out%estimate, 0.5_dp, 'binary accuracy')
   call near(out%se, 0.25_dp, 'accuracy standard error')
   out = loo_predictive_metric(xc, [0.0_dp, 0.0_dp, 1.0_dp, 1.0_dp], weights, 'balanced_acc')
   call near(out%estimate, 0.5_dp, 'balanced accuracy')
   call near(out%se, 0.176776695296637_dp, 'balanced accuracy standard error')

   out = crps(x(:, 1), x2(:, 1), y(1), permutations)
   call near(out%estimate, 0.3875_dp, 'single observation crps')
   call check(ieee_is_nan(out%se), 'single observation undefined standard error')
   out = scrps(x(:, 1), x2(:, 1), y(1), permutations)
   call near(out%estimate, -1.03703059242918_dp, 'single observation scrps')
   out = crps(x, x2, y)
   call check(out%status == 0, 'identity pairing accepted')
   out = scrps(0.0_dp * x, 0.0_dp * x2, y)
   call check(out%status == 2 .and. ieee_is_nan(out%estimate), 'degenerate scaled score')
   invalid = weights * 2.0_dp
   out = loo_predictive_metric(x, y, invalid)
   call check(out%status == 1, 'reject unnormalized weights')
   bad_perm = permutations
   bad_perm(1, 1) = bad_perm(2, 1)
   out = crps(x, x2, y, bad_perm)
   call check(out%status == 1, 'reject repeated permutation indices')
   out = loo_crps(x, x2, y, weights, joint(:, :, :1), permutations)
   call check(out%status == 1, 'reject missing joint permutation weights')
   out = loo_predictive_metric(xc, [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], weights, 'balanced_acc')
   call check(out%status == 2 .and. ieee_is_nan(out%estimate), 'missing class')
   out = loo_predictive_metric(x, y, weights, 'unknown')
   call check(out%status == 1, 'reject unknown metric')
   print *, 'loo metric reference and edge-case tests passed'
contains
   subroutine check(condition, label)
      logical, intent(in) :: condition !! Assertion that must evaluate to true.
      character(len = *), intent(in) :: label !! Diagnostic label for a failing assertion.
      if (.not. condition) then
         print *, 'FAIL: ', label
         error stop 1
      end if
   end subroutine check
   subroutine near(value, expected, label)
      real(dp), intent(in) :: value !! Computed scalar to compare with the R reference.
      real(dp), intent(in) :: expected !! Reference result from loo 2.10.1.
      character(len = *), intent(in) :: label !! Diagnostic label for a failing comparison.
      call check(abs(value - expected) < 2.0e-13_dp, label)
   end subroutine near
   subroutine near_vector(value, expected, label)
      real(dp), intent(in) :: value(:) !! Computed pointwise result.
      real(dp), intent(in) :: expected(:) !! R pointwise reference of the same length.
      character(len = *), intent(in) :: label !! Diagnostic label for a failing comparison.
      call check(all(abs(value - expected) < 2.0e-13_dp), label)
   end subroutine near_vector
end program test_metrics
