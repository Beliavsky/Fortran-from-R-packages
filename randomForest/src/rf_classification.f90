! SPDX-License-Identifier: GPL-2.0-or-later
module rf_classification
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use r_kinds, only : dp, i64
   use rf_rng, only : rf_rng_state, sample_indices, shuffle_real
   use rf_types, only : rf_options, rf_classification_forest
   use rf_tree_mod, only : build_class_tree, predict_class_tree, tree_var_used
   implicit none
   private

   public :: fit_classification, predict_classification
   public :: make_synthetic_class, fit_unsupervised

contains

   subroutine fit_classification(x, y, forest, ncat, options, class_weights, cutoff, case_weights, &
      strata, strata_sample_size, status, message)
      real(dp), intent(in) :: x(:,:)
      integer, intent(in) :: y(:)
      type(rf_classification_forest), intent(out) :: forest
      integer, intent(in), optional :: ncat(:)
      type(rf_options), intent(in), optional :: options
      real(dp), intent(in), optional :: class_weights(:), cutoff(:), case_weights(:)
      integer, intent(in), optional :: strata(:), strata_sample_size(:)
      integer, intent(out), optional :: status
      character(len=*), intent(out), optional :: message
      type(rf_options) :: opt
      type(rf_rng_state) :: rng
      integer, allocatable :: cats(:), sampled(:), inbag_count(:), obs_index(:), tree_pred(:), leaf(:)
      integer, allocatable :: prox_denom(:,:)
      real(dp), allocatable :: obs_weight(:), sample_weight(:), class_weight(:), prox_count(:,:)
      real(dp), allocatable :: imp_sum(:,:), imp_sumsq(:,:)
      logical, allocatable :: used(:)
      integer :: n, p, nclass, ntree, mtry, nodesize, maxnodes, sample_size
      integer :: t, i, k, n_unique, attempt, rc
      logical :: good_sample

      call set_ok(status, message)
      n = size(x, 1)
      p = size(x, 2)
      if (size(y) /= n .or. n <= 1 .or. p <= 0) then
         call set_error(1, 'x and y have incompatible or empty dimensions', status, message)
         return
      end if
      if (any(.not. ieee_is_finite(x))) then
         call set_error(2, 'x contains non-finite values, impute before fitting', status, message)
         return
      end if
      nclass = maxval(y)
      if (nclass < 2 .or. minval(y) < 1) then
         call set_error(3, 'classification labels must be consecutive positive integers starting at 1', status, message)
         return
      end if
      do i = 1, nclass
         if (count(y == i) == 0) then
            call set_error(4, 'classification labels must contain every class from 1 through max(y)', status, message)
            return
         end if
      end do

      call prepare_categories(x, ncat, cats, rc)
      if (rc /= 0) then
         call set_error(5, 'categorical predictors must be integer-coded from 1 through ncat', status, message)
         return
      end if
      opt = rf_options()
      if (present(options)) opt = options
      ntree = max(1, opt%ntree)
      mtry = opt%mtry
      if (mtry <= 0) mtry = max(1, int(sqrt(real(p, dp))))
      mtry = min(p, mtry)
      nodesize = opt%nodesize
      if (nodesize <= 0) nodesize = 1
      maxnodes = opt%maxnodes
      if (maxnodes <= 0) maxnodes = max(3, 2 * n + 1)
      sample_size = opt%sample_size
      if (sample_size <= 0) sample_size = n
      if (.not. opt%replace .and. sample_size > n .and. .not. present(strata)) then
         call set_error(6, 'sample_size cannot exceed n without replacement', status, message)
         return
      end if

      allocate(sample_weight(n))
      if (present(case_weights)) then
         if (size(case_weights) /= n .or. any(case_weights < 0.0_dp) .or. sum(case_weights) <= 0.0_dp) then
            call set_error(7, 'case_weights must be nonnegative, length n, and have positive sum', status, message)
            return
         end if
         sample_weight = case_weights / sum(case_weights)
      else
         sample_weight = 1.0_dp / real(n, dp)
      end if
      call make_class_weights(y, nclass, class_weights, class_weight, rc)
      if (rc /= 0) then
         call set_error(8, 'class_weights must be nonnegative, length nclass, and have positive sum', status, message)
         return
      end if

      if (present(strata) .neqv. present(strata_sample_size)) then
         call set_error(9, 'strata and strata_sample_size must be supplied together', status, message)
         return
      end if
      if (present(strata)) then
         if (size(strata) /= n .or. minval(strata) < 1 .or. maxval(strata) /= size(strata_sample_size)) then
            call set_error(10, 'strata must be length n and coded 1 through size(strata_sample_size)', status, message)
            return
         end if
         sample_size = sum(strata_sample_size)
      end if
      if (sample_size <= 0) then
         call set_error(11, 'the requested sample size must be positive', status, message)
         return
      end if

      forest%nclass = nclass
      forest%nvar = p
      forest%nobs = n
      allocate(forest%ncat(p), forest%trees(ntree), forest%cutoff(nclass))
      allocate(forest%oob_votes(n, nclass), forest%oob_count(n), forest%oob_prediction(n))
      allocate(forest%error_curve(ntree, nclass + 1), forest%importance_gini(p))
      forest%ncat = cats
      forest%oob_votes = 0
      forest%oob_count = 0
      forest%oob_prediction = 0
      forest%error_curve = 0.0_dp
      forest%importance_gini = 0.0_dp
      if (present(cutoff)) then
         if (size(cutoff) /= nclass .or. any(cutoff <= 0.0_dp)) then
            call set_error(12, 'cutoff must be positive and have length nclass', status, message)
            return
         end if
         forest%cutoff = cutoff
      else
         forest%cutoff = 1.0_dp / real(nclass, dp)
      end if
      if (opt%keep_inbag) then
         allocate(forest%inbag(n, ntree))
         forest%inbag = 0
      end if
      if (opt%importance) then
         allocate(forest%importance_accuracy(p, nclass + 1), forest%importance_sd(p, nclass + 1))
         forest%importance_accuracy = 0.0_dp
         forest%importance_sd = 0.0_dp
         allocate(imp_sum(p, nclass + 1), imp_sumsq(p, nclass + 1))
         imp_sum = 0.0_dp
         imp_sumsq = 0.0_dp
      end if
      if (opt%proximity) then
         allocate(prox_count(n, n), prox_denom(n, n))
         prox_count = 0.0_dp
         prox_denom = 0
      end if

      allocate(sampled(sample_size), inbag_count(n), obs_index(n), obs_weight(n))
      allocate(tree_pred(n), leaf(n), used(p))
      call rng%seed(opt%seed)

      do t = 1, ntree
         good_sample = .false.
         do attempt = 1, 30
            if (present(strata)) then
               call draw_stratified_sample(rng, strata, strata_sample_size, opt%replace, sample_weight, sampled, rc)
            else
               call sample_indices(rng, n, sample_size, opt%replace, sampled, sample_weight, rc)
            end if
            if (rc /= 0) then
               call set_error(13, 'sampling failed, check sampling weights and sample sizes', status, message)
               return
            end if
            inbag_count = 0
            do i = 1, sample_size
               inbag_count(sampled(i)) = inbag_count(sampled(i)) + 1
            end do
            good_sample = classes_in_sample(y, inbag_count, nclass) >= 2
            if (good_sample) exit
         end do
         if (.not. good_sample) then
            call set_error(14, 'fewer than two classes remained in the in-bag sample after 30 attempts', status, message)
            return
         end if
         if (opt%keep_inbag) forest%inbag(:, t) = inbag_count

         n_unique = count(inbag_count > 0)
         k = 0
         do i = 1, n
            if (inbag_count(i) > 0) then
               k = k + 1
               obs_index(k) = i
               obs_weight(k) = real(inbag_count(i), dp) * class_weight(y(i))
            end if
         end do
         call build_class_tree(x, y, obs_index(1:n_unique), obs_weight(1:n_unique), cats, nclass, &
            mtry, nodesize, maxnodes, rng, forest%trees(t))
         call predict_class_tree(forest%trees(t), x, cats, tree_pred, leaf)
         forest%importance_gini = forest%importance_gini + forest%trees(t)%impurity_decrease

         do i = 1, n
            if (inbag_count(i) == 0) then
               forest%oob_votes(i, tree_pred(i)) = forest%oob_votes(i, tree_pred(i)) + 1
               forest%oob_count(i) = forest%oob_count(i) + 1
            end if
         end do
         call update_oob_predictions(y, forest%oob_votes, forest%oob_count, forest%cutoff, rng, &
            forest%oob_prediction, forest%error_curve(t, :))

         if (opt%proximity) call accumulate_proximity(leaf, inbag_count, opt%oob_proximity, prox_count, prox_denom)

         if (opt%importance) then
            call tree_var_used(forest%trees(t), used)
            call classification_tree_importance(forest%trees(t), x, y, cats, inbag_count, tree_pred, used, &
               max(1, opt%n_perm), rng, imp_sum, imp_sumsq)
         end if
      end do

      forest%importance_gini = forest%importance_gini / real(ntree, dp)
      if (opt%importance) call finish_importance(imp_sum, imp_sumsq, ntree, forest%importance_accuracy, forest%importance_sd)
      if (opt%proximity) call finish_proximity(prox_count, prox_denom, ntree, opt%oob_proximity, forest%proximity)
   end subroutine fit_classification

   subroutine predict_classification(forest, x, prediction, probabilities, terminal_nodes)
      type(rf_classification_forest), intent(in) :: forest
      real(dp), intent(in) :: x(:,:)
      integer, intent(out) :: prediction(:)
      real(dp), intent(out), optional :: probabilities(:,:)
      integer, intent(out), optional :: terminal_nodes(:,:)
      integer, allocatable :: tree_pred(:), leaf(:), votes(:,:)
      integer :: n, t, i

      n = size(x, 1)
      if (size(prediction) /= n) error stop 'predict_classification: prediction has wrong size'
      if (size(x, 2) /= forest%nvar) error stop 'predict_classification: x has wrong number of variables'
      if (present(probabilities)) then
         if (size(probabilities, 1) /= n .or. size(probabilities, 2) /= forest%nclass) &
            error stop 'predict_classification: probabilities has wrong shape'
      end if
      if (present(terminal_nodes)) then
         if (size(terminal_nodes, 1) /= n .or. size(terminal_nodes, 2) /= size(forest%trees)) &
            error stop 'predict_classification: terminal_nodes has wrong shape'
      end if

      allocate(tree_pred(n), leaf(n), votes(n, forest%nclass))
      votes = 0
      do t = 1, size(forest%trees)
         call predict_class_tree(forest%trees(t), x, forest%ncat, tree_pred, leaf)
         do i = 1, n
            votes(i, tree_pred(i)) = votes(i, tree_pred(i)) + 1
         end do
         if (present(terminal_nodes)) terminal_nodes(:, t) = leaf
      end do
      do i = 1, n
         prediction(i) = cutoff_winner(votes(i, :), forest%cutoff)
      end do
      if (present(probabilities)) probabilities = real(votes, dp) / real(size(forest%trees), dp)
   end subroutine predict_classification


   subroutine fit_unsupervised(x, forest, ncat, options, status, message)
      real(dp), intent(in) :: x(:,:)
      type(rf_classification_forest), intent(out) :: forest
      integer, intent(in), optional :: ncat(:)
      type(rf_options), intent(in), optional :: options
      integer, intent(out), optional :: status
      character(len=*), intent(out), optional :: message
      type(rf_options) :: opt, one_opt
      type(rf_classification_forest) :: one
      real(dp), allocatable :: x2(:,:), synthetic(:,:), prox_count(:,:)
      integer, allocatable :: y2(:), cats(:), pred(:), leaf(:), denom(:,:), inbag_original(:)
      integer :: n, p, t, rc
      character(len=256) :: local_message

      call set_ok(status, message)
      n = size(x, 1)
      p = size(x, 2)
      if (n <= 1 .or. p <= 0) then
         call set_error(30, 'x must contain at least two observations and one variable', status, message)
         return
      end if
      call prepare_categories(x, ncat, cats, rc)
      if (rc /= 0) then
         call set_error(31, 'categorical predictors must be integer-coded from 1 through ncat', status, message)
         return
      end if
      opt = rf_options()
      if (present(options)) opt = options
      allocate(forest%ncat(p), forest%trees(max(1, opt%ntree)), forest%cutoff(2))
      allocate(forest%oob_votes(n, 2), forest%oob_count(n), forest%oob_prediction(n))
      allocate(forest%error_curve(max(1, opt%ntree), 3), forest%importance_gini(p))
      forest%nclass = 2
      forest%nvar = p
      forest%nobs = n
      forest%ncat = cats
      forest%cutoff = 0.5_dp
      forest%oob_votes = 0
      forest%oob_count = 0
      forest%oob_prediction = 0
      forest%error_curve = 0.0_dp
      forest%importance_gini = 0.0_dp
      if (opt%keep_inbag) then
         allocate(forest%inbag(n, max(1, opt%ntree)))
         forest%inbag = 0
      end if
      if (opt%importance) then
         allocate(forest%importance_accuracy(p, 3), forest%importance_sd(p, 3))
         forest%importance_accuracy = 0.0_dp
         forest%importance_sd = 0.0_dp
      end if
      if (opt%proximity) then
         allocate(prox_count(n, n), denom(n, n))
         prox_count = 0.0_dp
         denom = 0
      end if
      allocate(x2(2 * n, p), synthetic(n, p), y2(2 * n), pred(n), leaf(n), inbag_original(n))
      x2(1:n, :) = x
      y2(1:n) = 1
      y2(n + 1:2 * n) = 2

      do t = 1, size(forest%trees)
         call make_synthetic_class(x, synthetic, seed=int(modulo(opt%seed + int(t, i64) * 104729_i64, &
            int(huge(1), i64)), kind(1)))
         x2(n + 1:2 * n, :) = synthetic
         one_opt = opt
         one_opt%ntree = 1
         one_opt%seed = opt%seed + int(t, i64) * 65537_i64
         one_opt%proximity = .false.
         one_opt%keep_inbag = .true.
         call fit_classification(x2, y2, one, ncat=cats, options=one_opt, status=rc, message=local_message)
         if (rc /= 0) then
            call set_error(32, 'unsupervised tree fit failed: '//trim(local_message), status, message)
            return
         end if
         forest%trees(t) = one%trees(1)
         forest%importance_gini = forest%importance_gini + one%importance_gini
         if (opt%importance) then
            forest%importance_accuracy = forest%importance_accuracy + one%importance_accuracy
            forest%importance_sd = forest%importance_sd + one%importance_sd ** 2
         end if
         inbag_original = one%inbag(1:n, 1)
         if (opt%keep_inbag) forest%inbag(:, t) = inbag_original
         if (opt%proximity) then
            call predict_class_tree(forest%trees(t), x, cats, pred, leaf)
            call accumulate_proximity(leaf, inbag_original, opt%oob_proximity, prox_count, denom)
         end if
      end do
      forest%importance_gini = forest%importance_gini / real(size(forest%trees), dp)
      if (opt%importance) then
         forest%importance_accuracy = forest%importance_accuracy / real(size(forest%trees), dp)
         forest%importance_sd = sqrt(forest%importance_sd) / real(size(forest%trees), dp)
      end if
      if (opt%proximity) call finish_proximity(prox_count, denom, size(forest%trees), opt%oob_proximity, forest%proximity)
   end subroutine fit_unsupervised

   subroutine make_synthetic_class(x, synthetic, seed)
      real(dp), intent(in) :: x(:,:)
      real(dp), intent(out) :: synthetic(size(x, 1), size(x, 2))
      integer, intent(in), optional :: seed
      type(rf_rng_state) :: rng
      integer :: i, j

      if (present(seed)) then
         call rng%seed(int(seed, kind(rng%state)))
      else
         call rng%seed(int(1976, kind(rng%state)))
      end if
      do i = 1, size(x, 1)
         do j = 1, size(x, 2)
            synthetic(i, j) = x(rng%randint(size(x, 1)), j)
         end do
      end do
   end subroutine make_synthetic_class

   subroutine prepare_categories(x, ncat_in, cats, status)
      real(dp), intent(in) :: x(:,:)
      integer, intent(in), optional :: ncat_in(:)
      integer, allocatable, intent(out) :: cats(:)
      integer, intent(out) :: status
      integer :: j

      status = 0
      allocate(cats(size(x, 2)))
      if (present(ncat_in)) then
         if (size(ncat_in) /= size(x, 2) .or. any(ncat_in < 1) .or. any(ncat_in > 53)) then
            status = 1
            return
         end if
         cats = ncat_in
      else
         cats = 1
      end if
      do j = 1, size(cats)
         if (cats(j) > 1) then
            if (any(abs(x(:, j) - real(nint(x(:, j)), dp)) > 1.0e-10_dp)) then
               status = 2
               return
            end if
            if (minval(nint(x(:, j))) < 1 .or. maxval(nint(x(:, j))) > cats(j)) then
               status = 3
               return
            end if
         end if
      end do
   end subroutine prepare_categories

   subroutine make_class_weights(y, nclass, supplied, class_weight, status)
      integer, intent(in) :: y(:), nclass
      real(dp), intent(in), optional :: supplied(:)
      real(dp), allocatable, intent(out) :: class_weight(:)
      integer, intent(out) :: status
      integer :: c
      real(dp) :: s

      status = 0
      allocate(class_weight(nclass))
      if (present(supplied)) then
         if (size(supplied) /= nclass .or. any(supplied < 0.0_dp)) then
            status = 1
            return
         end if
         s = sum(supplied)
         if (s <= 0.0_dp) then
            status = 2
            return
         end if
         class_weight = supplied / s
      else
         do c = 1, nclass
            class_weight(c) = real(count(y == c), dp) / real(size(y), dp)
         end do
      end if
      do c = 1, nclass
         class_weight(c) = class_weight(c) * real(size(y), dp) / real(count(y == c), dp)
      end do
   end subroutine make_class_weights

   integer function classes_in_sample(y, inbag, nclass) result(n_present)
      integer, intent(in) :: y(:), inbag(:), nclass
      integer :: c

      n_present = 0
      do c = 1, nclass
         if (any(inbag > 0 .and. y == c)) n_present = n_present + 1
      end do
   end function classes_in_sample

   subroutine draw_stratified_sample(rng, strata, strata_sample_size, replace, weights, sampled, status)
      type(rf_rng_state), intent(inout) :: rng
      integer, intent(in) :: strata(:), strata_sample_size(:)
      logical, intent(in) :: replace
      real(dp), intent(in) :: weights(:)
      integer, intent(out) :: sampled(:)
      integer, intent(out) :: status
      integer, allocatable :: members(:), local_sample(:)
      real(dp), allocatable :: local_weights(:)
      integer :: s, nm, i, pos, rc

      status = 0
      pos = 1
      do s = 1, size(strata_sample_size)
         nm = count(strata == s)
         if (nm <= 0 .or. strata_sample_size(s) < 0) then
            status = 1
            return
         end if
         if (.not. replace .and. strata_sample_size(s) > nm) then
            status = 2
            return
         end if
         allocate(members(nm), local_weights(nm), local_sample(strata_sample_size(s)))
         members = pack([(i, i = 1, size(strata))], strata == s)
         local_weights = weights(members)
         call sample_indices(rng, nm, strata_sample_size(s), replace, local_sample, local_weights, rc)
         if (rc /= 0) then
            status = 3
            return
         end if
         do i = 1, strata_sample_size(s)
            sampled(pos) = members(local_sample(i))
            pos = pos + 1
         end do
         deallocate(members, local_weights, local_sample)
      end do
   end subroutine draw_stratified_sample

   subroutine update_oob_predictions(y, votes, counts, cutoff, rng, prediction, error)
      integer, intent(in) :: y(:), votes(:,:), counts(:)
      real(dp), intent(in) :: cutoff(:)
      type(rf_rng_state), intent(inout) :: rng
      integer, intent(out) :: prediction(:)
      real(dp), intent(out) :: error(:)
      integer :: i, c, n_used, nwrong, denom, wrong

      prediction = 0
      do i = 1, size(y)
         if (counts(i) > 0) prediction(i) = cutoff_winner_random(votes(i, :), cutoff, rng)
      end do
      n_used = count(counts > 0)
      if (n_used > 0) then
         nwrong = count(counts > 0 .and. prediction /= y)
         error(1) = real(nwrong, dp) / real(n_used, dp)
      else
         error(1) = 0.0_dp
      end if
      do c = 1, size(cutoff)
         denom = count(counts > 0 .and. y == c)
         if (denom > 0) then
            wrong = count(counts > 0 .and. y == c .and. prediction /= y)
            error(c + 1) = real(wrong, dp) / real(denom, dp)
         else
            error(c + 1) = 0.0_dp
         end if
      end do
   end subroutine update_oob_predictions

   integer function cutoff_winner(votes, cutoff) result(klass)
      integer, intent(in) :: votes(:)
      real(dp), intent(in) :: cutoff(:)
      integer :: c
      real(dp) :: score, best

      klass = 1
      best = real(votes(1), dp) / cutoff(1)
      do c = 2, size(votes)
         score = real(votes(c), dp) / cutoff(c)
         if (score > best) then
            best = score
            klass = c
         end if
      end do
   end function cutoff_winner

   integer function cutoff_winner_random(votes, cutoff, rng) result(klass)
      integer, intent(in) :: votes(:)
      real(dp), intent(in) :: cutoff(:)
      type(rf_rng_state), intent(inout) :: rng
      integer :: c, ntie
      real(dp) :: score, best, tol

      klass = 1
      best = real(votes(1), dp) / cutoff(1)
      ntie = 1
      do c = 2, size(votes)
         score = real(votes(c), dp) / cutoff(c)
         tol = 64.0_dp * epsilon(1.0_dp) * max(1.0_dp, abs(score), abs(best))
         if (score > best + tol) then
            best = score
            klass = c
            ntie = 1
         else if (abs(score - best) <= tol) then
            ntie = ntie + 1
            if (rng%uniform() < 1.0_dp / real(ntie, dp)) klass = c
         end if
      end do
   end function cutoff_winner_random

   subroutine accumulate_proximity(leaf, inbag, oob_only, count_matrix, denominator)
      integer, intent(in) :: leaf(:), inbag(:)
      logical, intent(in) :: oob_only
      real(dp), intent(inout) :: count_matrix(:,:)
      integer, intent(inout) :: denominator(:,:)
      integer :: i, j

      do i = 1, size(leaf)
         do j = i + 1, size(leaf)
            if (oob_only) then
               if (inbag(i) /= 0 .or. inbag(j) /= 0) cycle
               denominator(i, j) = denominator(i, j) + 1
               denominator(j, i) = denominator(i, j)
            end if
            if (leaf(i) == leaf(j)) then
               count_matrix(i, j) = count_matrix(i, j) + 1.0_dp
               count_matrix(j, i) = count_matrix(i, j)
            end if
         end do
      end do
   end subroutine accumulate_proximity

   subroutine finish_proximity(count_matrix, denominator, ntree, oob_only, proximity)
      real(dp), intent(in) :: count_matrix(:,:)
      integer, intent(in) :: denominator(:,:), ntree
      logical, intent(in) :: oob_only
      real(dp), allocatable, intent(out) :: proximity(:,:)
      integer :: i, j

      allocate(proximity(size(count_matrix, 1), size(count_matrix, 2)))
      proximity = 0.0_dp
      do i = 1, size(proximity, 1)
         proximity(i, i) = 1.0_dp
         do j = i + 1, size(proximity, 2)
            if (oob_only) then
               if (denominator(i, j) > 0) proximity(i, j) = count_matrix(i, j) / real(denominator(i, j), dp)
            else
               proximity(i, j) = count_matrix(i, j) / real(ntree, dp)
            end if
            proximity(j, i) = proximity(i, j)
         end do
      end do
   end subroutine finish_proximity

   subroutine classification_tree_importance(tree, x, y, ncat, inbag, base_pred, used, nperm, rng, imp_sum, imp_sumsq)
      use rf_types, only : rf_tree
      type(rf_tree), intent(in) :: tree
      real(dp), intent(in) :: x(:,:)
      integer, intent(in) :: y(:), ncat(:), inbag(:), base_pred(:), nperm
      logical, intent(in) :: used(:)
      type(rf_rng_state), intent(inout) :: rng
      real(dp), intent(inout) :: imp_sum(:,:), imp_sumsq(:,:)
      real(dp), allocatable :: xperm(:,:), vals(:), delta(:)
      integer, allocatable :: pred(:), oob_idx(:)
      integer :: m, k, i, c, noob, denom, base_correct, perm_correct
      real(dp) :: d

      noob = count(inbag == 0)
      if (noob <= 0) return
      allocate(oob_idx(noob), pred(size(y)), xperm(size(x, 1), size(x, 2)), vals(noob))
      allocate(delta(size(imp_sum, 2)))
      oob_idx = pack([(i, i = 1, size(y))], inbag == 0)

      do m = 1, size(used)
         if (.not. used(m)) cycle
         delta = 0.0_dp
         do k = 1, nperm
            xperm = x
            vals = x(oob_idx, m)
            call shuffle_real(rng, vals)
            xperm(oob_idx, m) = vals
            call predict_class_tree(tree, xperm, ncat, pred)
            base_correct = count(base_pred(oob_idx) == y(oob_idx))
            perm_correct = count(pred(oob_idx) == y(oob_idx))
            delta(1) = delta(1) + real(base_correct - perm_correct, dp) / real(noob, dp)
            do c = 1, size(delta) - 1
               denom = count(y(oob_idx) == c)
               if (denom > 0) then
                  base_correct = count(y(oob_idx) == c .and. base_pred(oob_idx) == y(oob_idx))
                  perm_correct = count(y(oob_idx) == c .and. pred(oob_idx) == y(oob_idx))
                  delta(c + 1) = delta(c + 1) + real(base_correct - perm_correct, dp) / real(denom, dp)
               end if
            end do
         end do
         delta = delta / real(nperm, dp)
         do c = 1, size(delta)
            d = delta(c)
            imp_sum(m, c) = imp_sum(m, c) + d
            imp_sumsq(m, c) = imp_sumsq(m, c) + d * d
         end do
      end do
   end subroutine classification_tree_importance

   subroutine finish_importance(sum_imp, sumsq_imp, ntree, mean_imp, sd_imp)
      real(dp), intent(in) :: sum_imp(:,:), sumsq_imp(:,:)
      integer, intent(in) :: ntree
      real(dp), intent(out) :: mean_imp(:,:), sd_imp(:,:)
      real(dp) :: v
      integer :: i, j

      mean_imp = sum_imp / real(ntree, dp)
      do i = 1, size(mean_imp, 1)
         do j = 1, size(mean_imp, 2)
            v = max(0.0_dp, sumsq_imp(i, j) / real(ntree, dp) - mean_imp(i, j) ** 2)
            sd_imp(i, j) = sqrt(v / real(ntree, dp))
         end do
      end do
   end subroutine finish_importance

   subroutine set_ok(status, message)
      integer, intent(out), optional :: status
      character(len=*), intent(out), optional :: message
      if (present(status)) status = 0
      if (present(message)) message = ''
   end subroutine set_ok

   subroutine set_error(code, text, status, message)
      integer, intent(in) :: code
      character(len=*), intent(in) :: text
      integer, intent(out), optional :: status
      character(len=*), intent(out), optional :: message
      if (present(status)) status = code
      if (present(message)) message = text
   end subroutine set_error

end module rf_classification
