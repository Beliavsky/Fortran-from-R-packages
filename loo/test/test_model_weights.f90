! SPDX-License-Identifier: GPL-3.0-or-later
! References generated using loo 2.10.1:
! stacking_weights(lpd,optim_control=list(reltol=1e-13,maxit=10000));
! pseudobma_weights(lpd,BB=FALSE); and set.seed(904) followed by
! pseudobma_weights(lpd,BB=TRUE,BB_n=4). Bootstrap rows are the exact draws from
! loo:::dirichlet_rng(4,rep(1,6)) under that seed, avoiding RNG compatibility assumptions.
program test_model_weights
   use r_kinds, only: dp
   use loo_model_weights
   use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_negative_inf, ieee_quiet_nan, ieee_is_nan
   implicit none
   real(dp) :: lpd(6, 3), weights(3), shifted_weights(3), bb(4, 6), objective, gap
   real(dp) :: two(2, 2), pair(2), sparse_bootstrap(2, 2), invalid(6, 3)
   integer :: status, iterations

   lpd = log(reshape([0.8_dp, 0.2_dp, 0.1_dp, 0.7_dp, 0.3_dp, 0.2_dp, &
      0.1_dp, 0.8_dp, 0.6_dp, 0.1_dp, 0.8_dp, 0.2_dp, &
      0.1_dp, 0.1_dp, 0.7_dp, 0.2_dp, 0.1_dp, 0.8_dp], [6, 3]))
   call stacking_weights(lpd, weights, status, objective, iterations, kkt_gap = gap)
   call check(status == 0 .and. gap <= 1.0e-10_dp, 'stacking convergence certificate')
   call check(iterations > 0, 'stacking iteration count')
   call near_vector(weights, [0.415385134208523_dp, 0.480506922832286_dp, 0.104107942959192_dp], &
      3.0e-6_dp, 'stacking R reference')
   call check(abs(sum(weights) - 1.0_dp) < 2.0e-15_dp .and. minval(weights) >= 0.0_dp, 'simplex constraint')
   call check(abs(objective - sum(log(matmul(exp(lpd), weights)))) < 2.0e-14_dp, 'stacking objective')
   call stacking_weights(lpd + 1000.0_dp, shifted_weights, status)
   call check(status == 0, 'large log density convergence')
   call near_vector(weights, shifted_weights, 3.0e-10_dp, 'log-scale shift invariance')
   call model_weights(lpd, shifted_weights, status = status)
   call near_vector(weights, shifted_weights, 2.0e-14_dp, 'stacking wrapper')
   call stacking_weights(lpd, weights, status, max_iterations = 1)
   call check(status == 2, 'report nonconvergence')

   call pseudobma_weights(lpd, weights, status = status)
   call check(status == 0, 'pseudo-BMA status')
   call near_vector(weights, [0.432989690721650_dp, 0.494845360824742_dp, 0.0721649484536082_dp], &
      3.0e-14_dp, 'pseudo-BMA reference')
   bb = reshape([0.205040841908457_dp, 0.129604225518116_dp, 0.270265992691146_dp, &
      0.0648563886392211_dp, 0.0106475941993484_dp, 0.380094682496800_dp, &
      0.0619939397708617_dp, 0.127972629112624_dp, 0.363965606510990_dp, &
      0.267779386590858_dp, 0.191559032563588_dp, 0.0717630442329597_dp, &
      0.0495946292630304_dp, 0.0550137226156912_dp, 0.181294041414920_dp, &
      0.678395102491341_dp, 0.197775297381351_dp, 0.0809523966663572_dp, &
      0.217387309697147_dp, 0.0548897539635107_dp, 0.172976030736823_dp, &
      0.0865555861121779_dp, 0.0774996838623377_dp, 0.00212308156034411_dp], [4, 6])
   call pseudobma_weights(lpd, weights, bb, status)
   call check(status == 0, 'pseudo-BMA+ status')
   call near_vector(weights, [0.476406657124527_dp, 0.433858017576657_dp, 0.0897353252988165_dp], &
      3.0e-14_dp, 'pseudo-BMA+ reference')
   call model_weights(lpd, shifted_weights, 'pseudobma+', bb, status)
   call near_vector(weights, shifted_weights, 2.0e-14_dp, 'pseudo-BMA+ wrapper')
   call model_weights(lpd, weights, 'pseudobma+', status = status)
   call check(status == 1 .and. all(ieee_is_nan(weights)), 'require bootstrap for plus')
   call pseudobma_weights(lpd, weights, bb * 2.0_dp, status)
   call check(status == 1, 'reject unnormalized bootstrap')

   two = log(reshape([0.9_dp, 0.8_dp, 0.1_dp, 0.2_dp], [2, 2]))
   call stacking_weights(two, pair, status, kkt_gap = gap)
   call check(status == 0, 'boundary optimum convergence')
   call near_vector(pair, [1.0_dp, 0.0_dp], 1.0e-15_dp, 'boundary optimum')
   two = -5000.0_dp
   call stacking_weights(two, pair, status)
   call check(status == 0, 'identical model convergence')
   call near_vector(pair, [0.5_dp, 0.5_dp], 1.0e-15_dp, 'identical models')
   two = ieee_value(0.0_dp, ieee_negative_inf)
   two(1, 1) = 0.0_dp
   two(2, 2) = 0.0_dp
   call stacking_weights(two, pair, status)
   call check(status == 0, 'disjoint support')
   call near_vector(pair, [0.5_dp, 0.5_dp], 1.0e-15_dp, 'disjoint support optimum')
   call pseudobma_weights(two, pair, status = status)
   call check(status == 1 .and. all(ieee_is_nan(pair)), 'all totals negative infinity')
   sparse_bootstrap = reshape([1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [2, 2])
   call pseudobma_weights(two, pair, sparse_bootstrap, status)
   call check(status == 0, 'sparse bootstrap handles zero times infinity')
   call near_vector(pair, [0.5_dp, 0.5_dp], 1.0e-15_dp, 'sparse bootstrap reference')
   invalid = lpd
   invalid(1, :) = ieee_value(0.0_dp, ieee_negative_inf)
   call stacking_weights(invalid, weights, status)
   call check(status == 1, 'reject all-infinite observation')
   invalid = lpd
   invalid(2, 2) = ieee_value(0.0_dp, ieee_quiet_nan)
   call stacking_weights(invalid, weights, status)
   call check(status == 1, 'reject NaN')
   print *, 'loo model-weight reference and edge-case tests passed'
contains
   subroutine check(condition, label)
      logical, intent(in) :: condition !! Assertion that must evaluate to true.
      character(len = *), intent(in) :: label !! Diagnostic label for a failing assertion.
      if (.not. condition) then
         print *, 'FAIL: ', label
         error stop 1
      end if
   end subroutine check
   subroutine near_vector(value, expected, tolerance, label)
      real(dp), intent(in) :: value(:) !! Computed model-weight vector.
      real(dp), intent(in) :: expected(:) !! Reference model-weight vector of the same length.
      real(dp), intent(in) :: tolerance !! Maximum absolute componentwise discrepancy.
      character(len = *), intent(in) :: label !! Diagnostic label for a failing comparison.
      call check(all(abs(value - expected) < tolerance), label)
   end subroutine near_vector
end program test_model_weights
