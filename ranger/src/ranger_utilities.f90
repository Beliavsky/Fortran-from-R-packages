! SPDX-License-Identifier: GPL-3.0-only
! Modern Fortran translation of ranger 0.18.0 computational code.
! Upstream C++ core: Copyright (c) 2014-2018 Marvin N. Wright; MIT licensed.
! The upstream R package is GPL-3. See NOTICE.md for full provenance.
module ranger_utilities
   use r_kinds, only : dp, i64
   use ranger_ij_calibration, only : calibrate_ij_variances
   use ranger_types, only : ranger_regression_forest, ranger_probability_forest, ranger_tree, RANGER_INTERIOR
   implicit none
   private

   public :: importance_pvalues_janitza, importance_pvalues_altmann
   public :: infinitesimal_jackknife, hierarchical_shrink_regression, hierarchical_shrink_probability
   public :: case_specific_weights, tree_sizes, variable_usage

contains

   subroutine importance_pvalues_janitza(importance, pvalue, status)
      real(dp), intent(in) :: importance(:)
      real(dp), intent(out) :: pvalue(size(importance))
      integer, intent(out), optional :: status
      real(dp), allocatable :: null_values(:)
      integer :: nneg, nzero, nnull, i, j, nle

      if (present(status)) status = 0
      nneg = count(importance < 0.0_dp)
      nzero = count(abs(importance) <= 0.0_dp)
      if (nneg == 0) then
         pvalue = 1.0_dp
         if (present(status)) status = 1
         return
      end if
      nnull = 2 * nneg + nzero
      allocate(null_values(nnull))
      j = 0
      do i = 1, size(importance)
         if (importance(i) < 0.0_dp) then
            j = j + 1
            null_values(j) = importance(i)
            j = j + 1
            null_values(j) = -importance(i)
         else if (abs(importance(i)) <= 0.0_dp) then
            j = j + 1
            null_values(j) = 0.0_dp
         end if
      end do
      do i = 1, size(importance)
         nle = count(null_values < importance(i))
         pvalue(i) = 1.0_dp - real(nle, dp) / real(nnull, dp)
         pvalue(i) = max(0.0_dp, min(1.0_dp, pvalue(i)))
      end do
   end subroutine importance_pvalues_janitza

   subroutine importance_pvalues_altmann(observed, null_importance, pvalue)
      real(dp), intent(in) :: observed(:), null_importance(:,:)
      real(dp), intent(out) :: pvalue(size(observed))
      integer :: i, nperm

      if (size(null_importance, 1) /= size(observed)) &
         error stop 'importance_pvalues_altmann: first dimension must match observed importance'
      nperm = size(null_importance, 2)
      if (nperm <= 0) error stop 'importance_pvalues_altmann: at least one permutation is required'
      do i = 1, size(observed)
         pvalue(i) = real(count(null_importance(i, :) >= observed(i)) + 1, dp) / real(nperm + 1, dp)
      end do
   end subroutine importance_pvalues_altmann

   subroutine infinitesimal_jackknife(pred, inbag, mean_prediction, variance, calibrate, seed)
      real(dp), intent(in) :: pred(:,:)
      integer, intent(in) :: inbag(:,:)
      real(dp), intent(out) :: mean_prediction(size(pred, 1)), variance(size(pred, 1))
      logical, intent(in), optional :: calibrate
      integer(i64), intent(in), optional :: seed
      real(dp), allocatable :: centered(:,:), navg(:), covariance(:,:)
      real(dp) :: nvar, bootvar, inflation, sample_fraction
      real(dp), allocatable :: variance_calibrated(:)
      integer :: ntrain, ntest, ntree, i, j, b
      logical :: do_calibrate, no_replacement

      ntest = size(pred, 1)
      ntree = size(pred, 2)
      ntrain = size(inbag, 1)
      if (size(inbag, 2) /= ntree .or. ntree <= 1) &
         error stop 'infinitesimal_jackknife: pred and inbag have incompatible dimensions'
      allocate(centered(ntest, ntree), navg(ntrain), covariance(ntrain, ntest))
      mean_prediction = sum(pred, dim=2) / real(ntree, dp)
      do b = 1, ntree
         centered(:, b) = pred(:, b) - mean_prediction
      end do
      navg = sum(real(inbag, dp), dim=2) / real(ntree, dp)
      covariance = 0.0_dp
      do i = 1, ntrain
         do j = 1, ntest
            do b = 1, ntree
               covariance(i, j) = covariance(i, j) + &
                  (real(inbag(i, b), dp) - navg(i)) * centered(j, b)
            end do
            covariance(i, j) = covariance(i, j) / real(ntree, dp)
         end do
      end do
      variance = sum(covariance * covariance, dim=1)

      nvar = 0.0_dp
      do i = 1, ntrain
         nvar = nvar + sum(real(inbag(i, :), dp) ** 2) / real(ntree, dp) - navg(i) ** 2
      end do
      nvar = nvar / real(ntrain, dp)
      do j = 1, ntest
         bootvar = sum(centered(j, :) ** 2) / real(ntree, dp)
         variance(j) = variance(j) - real(ntrain, dp) * nvar * bootvar / real(ntree, dp)
      end do

      no_replacement = maxval(inbag) <= 1
      if (no_replacement) then
         sample_fraction = sum(real(inbag, dp)) / real(ntrain * ntree, dp)
         if (sample_fraction < 1.0_dp) then
            inflation = 1.0_dp / (1.0_dp - sample_fraction) ** 2
            variance = inflation * variance
         end if
      end if
      do_calibrate = .false.
      if (present(calibrate)) do_calibrate = calibrate
      if (do_calibrate .and. ntest > 20) then
         allocate(variance_calibrated(ntest))
         if (present(seed)) then
            call calibrate_ij_variances(pred, inbag, variance, variance_calibrated, seed=seed)
         else
            call calibrate_ij_variances(pred, inbag, variance, variance_calibrated)
         end if
         variance = variance_calibrated
      end if
   end subroutine infinitesimal_jackknife

   subroutine hierarchical_shrink_regression(forest, lambda)
      type(ranger_regression_forest), intent(inout) :: forest
      real(dp), intent(in) :: lambda
      integer :: t
      if (lambda < 0.0_dp) error stop 'hierarchical_shrink_regression: lambda must be nonnegative'
      do t = 1, size(forest%trees)
         call shrink_reg_tree(forest%trees(t), 1, 0, 0.0_dp, 0.0_dp, lambda)
      end do
   end subroutine hierarchical_shrink_regression

   recursive subroutine shrink_reg_tree(tree, node, parent_n, parent_prediction, cumulative, lambda)
      type(ranger_tree), intent(inout) :: tree
      integer, intent(in) :: node, parent_n
      real(dp), intent(in) :: parent_prediction, cumulative, lambda
      real(dp) :: current

      if (node == 1) then
         current = tree%node_mean(node)
      else
         current = cumulative + (tree%node_mean(node) - parent_prediction) / &
            (1.0_dp + lambda / real(max(1, parent_n), dp))
      end if
      if (tree%status(node) /= RANGER_INTERIOR) then
         tree%node_mean(node) = current
      else
         call shrink_reg_tree(tree, tree%left(node), tree%node_n(node), tree%node_mean(node), current, lambda)
         call shrink_reg_tree(tree, tree%right(node), tree%node_n(node), tree%node_mean(node), current, lambda)
      end if
   end subroutine shrink_reg_tree

   subroutine hierarchical_shrink_probability(forest, lambda)
      type(ranger_probability_forest), intent(inout) :: forest
      real(dp), intent(in) :: lambda
      real(dp), allocatable :: parent(:), cumulative(:)
      integer :: t
      if (lambda < 0.0_dp) error stop 'hierarchical_shrink_probability: lambda must be nonnegative'
      allocate(parent(forest%nclass), cumulative(forest%nclass))
      do t = 1, size(forest%trees)
         parent = 0.0_dp
         cumulative = 0.0_dp
         call shrink_prob_tree(forest%trees(t), 1, 0, parent, cumulative, lambda)
      end do
   end subroutine hierarchical_shrink_probability

   recursive subroutine shrink_prob_tree(tree, node, parent_n, parent_prediction, cumulative, lambda)
      type(ranger_tree), intent(inout) :: tree
      integer, intent(in) :: node, parent_n
      real(dp), intent(in) :: parent_prediction(:), cumulative(:), lambda
      real(dp), allocatable :: current(:)

      allocate(current(size(cumulative)))
      if (node == 1) then
         current = tree%class_prob(:, node)
      else
         current = cumulative + (tree%class_prob(:, node) - parent_prediction) / &
            (1.0_dp + lambda / real(max(1, parent_n), dp))
      end if
      if (tree%status(node) /= RANGER_INTERIOR) then
         tree%class_prob(:, node) = max(current, 0.0_dp)
         if (sum(tree%class_prob(:, node)) > 0.0_dp) &
            tree%class_prob(:, node) = tree%class_prob(:, node) / sum(tree%class_prob(:, node))
      else
         call shrink_prob_tree(tree, tree%left(node), tree%node_n(node), tree%class_prob(:, node), current, lambda)
         call shrink_prob_tree(tree, tree%right(node), tree%node_n(node), tree%class_prob(:, node), current, lambda)
      end if
   end subroutine shrink_prob_tree

   subroutine case_specific_weights(training_terminal, test_terminal, weights)
      integer, intent(in) :: training_terminal(:,:), test_terminal(:)
      real(dp), intent(out) :: weights(size(training_terminal, 1))
      integer :: i, ntree
      real(dp) :: total

      ntree = size(training_terminal, 2)
      if (size(test_terminal) /= ntree) error stop 'case_specific_weights: terminal-node tree counts differ'
      do i = 1, size(training_terminal, 1)
         weights(i) = real(count(training_terminal(i, :) == test_terminal), dp)
      end do
      total = sum(weights)
      if (total > 0.0_dp) then
         weights = weights / total
      else
         weights = 1.0_dp / real(size(weights), dp)
      end if
   end subroutine case_specific_weights

   subroutine tree_sizes(trees, sizes)
      type(ranger_tree), intent(in) :: trees(:)
      integer, intent(out) :: sizes(size(trees))
      integer :: t
      do t = 1, size(trees)
         sizes(t) = count(trees(t)%status(1:trees(t)%n_nodes) /= RANGER_INTERIOR)
      end do
   end subroutine tree_sizes

   subroutine variable_usage(trees, nvar, counts)
      type(ranger_tree), intent(in) :: trees(:)
      integer, intent(in) :: nvar
      integer, intent(out) :: counts(nvar)
      integer :: t, node, v
      counts = 0
      do t = 1, size(trees)
         do node = 1, trees(t)%n_nodes
            if (trees(t)%status(node) /= RANGER_INTERIOR) cycle
            v = trees(t)%split_var(node)
            if (v >= 1 .and. v <= nvar) counts(v) = counts(v) + 1
         end do
      end do
   end subroutine variable_usage

end module ranger_utilities
