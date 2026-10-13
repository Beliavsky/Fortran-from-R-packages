! SPDX-License-Identifier: GPL-3.0-only
program test_properties
   use prroc_api
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_value, ieee_quiet_nan, ieee_positive_inf
   implicit none
   type(prroc_result) :: result, baseline, other
   real(dp) :: wp(3), wn(3), sp(3), sn(3), factor, pair_area, nan, inf
   integer :: status, i, j, checks

   checks = 0
   nan = ieee_value(0.0_dp, ieee_quiet_nan)
   inf = ieee_value(0.0_dp, ieee_positive_inf)
   call pr_curve([0.0_dp, 0.0_dp], result, scores_class1 = [1.0_dp, 1.0_dp], status = status)
   call require(status == prroc_ok, 'reverse PR status')
   call close(result%auc, 1.0_dp - log(2.0_dp), 'nonlinear reverse PR')
   call close(result%auc_davis_goadrich, 7.0_dp / 24.0_dp, 'Davis-Goadrich reverse PR')
   call require(abs(result%auc - 0.25_dp) > 0.05_dp, 'PR is not threshold trapezoids')
   call require(size(result%curve, 1) == 0, 'curve omitted by default')

   sp = [5.0_dp, 3.0_dp, 1.0_dp]
   sn = [4.0_dp, 3.0_dp, 0.0_dp]
   wp = [3.0_dp, 2.0_dp, 1.0_dp]
   wn = [2.0_dp, 5.0_dp, 4.0_dp]
   pair_area = 0.0_dp
   do i = 1, 3
      do j = 1, 3
         if (sp(i) > sn(j)) pair_area = pair_area + wp(i) * wn(j)
         if (sp(i) == sn(j)) pair_area = pair_area + 0.5_dp * wp(i) * wn(j)
      end do
   end do
   pair_area = pair_area / sum(wp) / sum(wn)
   call roc_curve(sp, result, sn, wp, wn, curve = .true., status = status)
   call require(status == prroc_ok, 'weighted ROC status')
   call close(result%auc, pair_area, 'independent weighted pairwise ROC identity')
   call roc_curve(-sp, other, -sn, wp, wn)
   call close(other%auc + result%auc, 1.0_dp, 'score reversal ROC identity')
   call pr_curve(sp, baseline, sn, wp, wn, max_compute = .true., min_compute = .true., rand_compute = .true.)
   do i = 1, 2
      factor = 1.0e250_dp
      if (i == 2) factor = 1.0e-250_dp
      call pr_curve(sp, result, sn, wp * factor, wn * factor, max_compute = .true., &
         min_compute = .true., rand_compute = .true., status = status)
      call require(status == prroc_ok, 'extreme common PR scale status')
      call close(result%auc, baseline%auc, 'PR scale invariance')
      call close(result%maximum%auc, baseline%maximum%auc, 'maximum scale invariance')
      call close(result%minimum%auc, baseline%minimum%auc, 'minimum scale invariance')
      call close(result%random%auc, baseline%random%auc, 'prevalence scale invariance')
      call roc_curve(sp, result, sn, wp * factor, wn * factor, status = status)
      call require(status == prroc_ok, 'extreme common ROC scale status')
      call close(result%auc, pair_area, 'ROC scale invariance')
   end do

   call pr_curve([0.0_dp], result, [1.0_dp], [1.0e-20_dp], [1.0_dp], status = status)
   call require(status == prroc_ok, 'tiny positive class status')
   call require(abs(result%auc / 5.0e-21_dp - 1.0_dp) < 2.0e-14_dp, 'small increment avoids cancellation')
   call pr_curve([2.0_dp], result, [1.0_dp], [1.0e-250_dp], [1.0e250_dp], status = status)
   call require(status == prroc_ok, '500-order imbalance PR status')
   call close(result%auc, 1.0_dp, '500-order imbalance perfect PR')
   call roc_curve([2.0_dp], result, [1.0_dp], [1.0e-250_dp], [1.0e250_dp], status = status)
   call require(status == prroc_ok, '500-order imbalance ROC status')
   call close(result%auc, 1.0_dp, '500-order imbalance perfect ROC')
   call pr_curve([0.0_dp], result, [0.0_dp], [1.0e308_dp], [1.0e308_dp], rand_compute = .true.)
   call close(result%auc, 0.5_dp, 'overflow-safe tied PR')
   call close(result%random%auc, 0.5_dp, 'overflow-safe prevalence')

   call pr_curve([3.0_dp, 1.0_dp], result, [4.0_dp, 2.0_dp, 0.0_dp], [.25_dp, .8_dp], &
      [.2_dp, .6_dp, .4_dp], rand_compute = .true.)
   call close(result%random%auc, 1.05_dp / 2.25_dp, 'unequal lengths use actual prevalence')
   call pr_curve([9.0_dp, 2.0_dp, -9.0_dp], result, [8.0_dp, 1.0_dp, -8.0_dp], &
      [0.0_dp, 1.0_dp, 0.0_dp], [0.0_dp, 1.0_dp, 0.0_dp], curve = .true.)
   call pr_curve([2.0_dp], other, [1.0_dp], [1.0_dp], [1.0_dp], curve = .true.)
   call close(result%auc, other%auc, 'zero-weight rows do not alter PR')
   call require(all(shape(result%curve) == shape(other%curve)), 'zero-weight thresholds omitted')
   call require(all(result%curve == other%curve), 'zero-weight curves agree')
   call pr_curve([2.0_dp], result, [1.0_dp], [1.0_dp], [0.0_dp], status = status)
   call require(status == prroc_ok, 'positive-only PR status')
   call close(result%auc, 1.0_dp, 'positive-only PR')
   call roc_curve([2.0_dp], result, [1.0_dp], [1.0_dp], [0.0_dp], status = status)
   call require(status == prroc_degenerate .and. ieee_is_nan(result%auc), 'one-class ROC rejected')
   call pr_curve([2.0_dp], result, weights_class0 = [0.0_dp], status = status)
   call require(status == prroc_degenerate .and. ieee_is_nan(result%auc), 'no-positive PR rejected')
   call pr_curve([2.0_dp], result, [1.0_dp], [0.0_dp], [0.0_dp], status = status)
   call require(status == prroc_degenerate, 'all-zero masses rejected')
   call pr_curve([real(dp) ::], result, [real(dp) ::], status = status)
   call require(status == prroc_invalid_input, 'empty input rejected')
   call pr_curve([2.0_dp], result, [real(dp) ::], status = status)
   call require(status == prroc_ok, 'empty negative vector permits PR')
   call close(result%auc, 1.0_dp, 'empty negative class PR')

   call pr_curve([nan], result, status = status)
   call require(status == prroc_invalid_input, 'NaN score rejected')
   call pr_curve([2.0_dp], result, weights_class0 = [nan], status = status)
   call require(status == prroc_invalid_input, 'NaN weight rejected')
   call pr_curve([2.0_dp], result, weights_class0 = [inf], status = status)
   call require(status == prroc_invalid_input, 'infinite weight rejected')
   call pr_curve([2.0_dp], result, weights_class0 = [-1.0_dp], status = status)
   call require(status == prroc_invalid_input, 'negative weight rejected')
   call pr_curve([2.0_dp], result, weights_class0 = [2.0_dp], status = status)
   call require(status == prroc_invalid_input, 'negative default complement rejected')
   call pr_curve([2.0_dp], result, weights_class0 = [1.0_dp, 1.0_dp], status = status)
   call require(status == prroc_invalid_input, 'weight shape validated')
   call pr_curve([2.0_dp], result, [1.0_dp, 0.0_dp], [.5_dp], status = status)
   call require(status == prroc_invalid_input, 'complement shape validated')
   call pr_curve([2.0_dp, 1.0_dp], result, sorted = .true., status = status)
   call require(status == prroc_invalid_input, 'sorted precondition validated')
   call pr_curve([2.0_dp], result, min_step_size = 0.0_dp, status = status)
   call require(status == prroc_invalid_input, 'zero curve step rejected')
   call pr_curve([2.0_dp], result, min_step_size = nan, status = status)
   call require(status == prroc_invalid_input, 'NaN curve step rejected')
   call pr_curve([2.0_dp], result, max_curve_points = 1, status = status)
   call require(status == prroc_invalid_input, 'invalid curve limit rejected')
   call pr_curve([inf], result, [-inf], status = status)
   call require(status == prroc_ok, 'infinite thresholds supported')
   call close(result%auc, 1.0_dp, 'infinite threshold PR')
   call roc_curve([inf], result, [-inf], status = status)
   call close(result%auc, 1.0_dp, 'infinite threshold ROC')
   call pr_curve([0.0_dp], result, [1.0_dp], curve = .true., min_step_size = 1.0e-20_dp, &
      max_curve_points = 10, status = status)
   call require(status == prroc_curve_limit, 'curve allocation bounded')
   call close(result%auc, 1.0_dp - log(2.0_dp), 'curve limit preserves AUC')
   call require(size(result%curve, 1) == 0, 'oversize curve is empty')
   call pr_curve([0.0_dp], result, [1.0_dp], dg_compute = .false.)
   call require(ieee_is_nan(result%auc_davis_goadrich), 'DG disabled')
   print *, 'PRROC independent property assertions:', checks
contains
   subroutine require(condition, description)
      logical, intent(in) :: condition !! Assertion that must hold for the named independent property.
      character(*), intent(in) :: description !! Diagnostic label printed if the assertion fails.
      if (.not. condition) then
         print *, description
         error stop 'PRROC property failure'
      end if
      checks = checks + 1
   end subroutine require

   subroutine close(actual, expected, description)
      real(dp), intent(in) :: actual !! Computed scalar required to be finite and close to the target.
      real(dp), intent(in) :: expected !! Independently derived identity or baseline scalar.
      character(*), intent(in) :: description !! Mathematical or invariance property being checked.
      call require(.not. ieee_is_nan(actual), description)
      call require(abs(actual - expected) <= 3.0e-13_dp, description)
   end subroutine close
end program test_properties
