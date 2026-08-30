! SPDX-License-Identifier: GPL-2.0-or-later
module rf_regression
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use r_kinds, only : dp
   use rf_rng, only : rf_rng_state, sample_indices, shuffle_real
   use rf_types, only : rf_options, rf_regression_forest, rf_tree
   use rf_tree_mod, only : build_reg_tree, predict_reg_tree, tree_var_used
   implicit none
   private

   public :: fit_regression, predict_regression

contains

   subroutine fit_regression(x, y, forest, ncat, options, case_weights, status, message)
      real(dp), intent(in) :: x(:,:), y(:)
      type(rf_regression_forest), intent(out) :: forest
      integer, intent(in), optional :: ncat(:)
      type(rf_options), intent(in), optional :: options
      real(dp), intent(in), optional :: case_weights(:)
      integer, intent(out), optional :: status
      character(len=*), intent(out), optional :: message
      type(rf_options) :: opt
      type(rf_rng_state) :: rng
      integer, allocatable :: cats(:), sampled(:), inbag_count(:), obs_index(:), leaf(:), prox_denom(:,:)
      real(dp), allocatable :: obs_weight(:), sample_weight(:), tree_pred(:), prox_count(:,:)
      real(dp), allocatable :: imp_sum(:), imp_sumsq(:)
      logical, allocatable :: used(:)
      integer :: n, p, ntree, mtry, nodesize, maxnodes, sample_size
      integer :: t, i, k, n_unique, rc, n_oob
      real(dp) :: mse, intercept, slope

      call set_ok(status, message)
      n = size(x, 1)
      p = size(x, 2)
      if (size(y) /= n .or. n <= 1 .or. p <= 0) then
         call set_error(1, 'x and y have incompatible or empty dimensions', status, message)
         return
      end if
      if (any(.not. ieee_is_finite(x)) .or. any(.not. ieee_is_finite(y))) then
         call set_error(2, 'x and y must contain only finite values', status, message)
         return
      end if
      call prepare_categories(x, ncat, cats, rc)
      if (rc /= 0) then
         call set_error(3, 'categorical predictors must be integer-coded from 1 through ncat', status, message)
         return
      end if

      opt = rf_options()
      if (present(options)) opt = options
      ntree = max(1, opt%ntree)
      mtry = opt%mtry
      if (mtry <= 0) mtry = max(1, p / 3)
      mtry = min(p, mtry)
      nodesize = opt%nodesize
      if (nodesize <= 0) nodesize = 5
      maxnodes = opt%maxnodes
      if (maxnodes <= 0) maxnodes = max(3, 2 * n + 1)
      sample_size = opt%sample_size
      if (sample_size <= 0) sample_size = n
      if (.not. opt%replace .and. sample_size > n) then
         call set_error(4, 'sample_size cannot exceed n without replacement', status, message)
         return
      end if

      allocate(sample_weight(n))
      if (present(case_weights)) then
         if (size(case_weights) /= n .or. any(case_weights < 0.0_dp) .or. sum(case_weights) <= 0.0_dp) then
            call set_error(5, 'case_weights must be nonnegative, length n, and have positive sum', status, message)
            return
         end if
         sample_weight = case_weights / sum(case_weights)
      else
         sample_weight = 1.0_dp / real(n, dp)
      end if

      forest%nvar = p
      forest%nobs = n
      allocate(forest%ncat(p), forest%trees(ntree), forest%oob_prediction(n), forest%oob_count(n))
      allocate(forest%mse_curve(ntree), forest%importance_gini(p))
      forest%ncat = cats
      forest%oob_prediction = 0.0_dp
      forest%oob_count = 0
      forest%mse_curve = 0.0_dp
      forest%importance_gini = 0.0_dp
      if (opt%keep_inbag) then
         allocate(forest%inbag(n, ntree))
         forest%inbag = 0
      end if
      if (opt%importance) then
         allocate(forest%importance_accuracy(p), forest%importance_sd(p), imp_sum(p), imp_sumsq(p))
         forest%importance_accuracy = 0.0_dp
         forest%importance_sd = 0.0_dp
         imp_sum = 0.0_dp
         imp_sumsq = 0.0_dp
      end if
      if (opt%proximity) then
         allocate(prox_count(n, n), prox_denom(n, n))
         prox_count = 0.0_dp
         prox_denom = 0
      end if

      allocate(sampled(sample_size), inbag_count(n), obs_index(n), obs_weight(n), tree_pred(n), leaf(n), used(p))
      call rng%seed(opt%seed)
      forest%bias_intercept = 0.0_dp
      forest%bias_slope = 1.0_dp

      do t = 1, ntree
         call sample_indices(rng, n, sample_size, opt%replace, sampled, sample_weight, rc)
         if (rc /= 0) then
            call set_error(6, 'sampling failed, check sampling weights and sample size', status, message)
            return
         end if
         inbag_count = 0
         do i = 1, sample_size
            inbag_count(sampled(i)) = inbag_count(sampled(i)) + 1
         end do
         if (opt%keep_inbag) forest%inbag(:, t) = inbag_count
         n_unique = count(inbag_count > 0)
         k = 0
         do i = 1, n
            if (inbag_count(i) > 0) then
               k = k + 1
               obs_index(k) = i
               obs_weight(k) = real(inbag_count(i), dp)
            end if
         end do

         call build_reg_tree(x, y, obs_index(1:n_unique), obs_weight(1:n_unique), cats, mtry, nodesize, maxnodes, &
            rng, forest%trees(t))
         call predict_reg_tree(forest%trees(t), x, cats, tree_pred, leaf)
         forest%importance_gini = forest%importance_gini + forest%trees(t)%impurity_decrease

         do i = 1, n
            if (inbag_count(i) == 0) then
               forest%oob_count(i) = forest%oob_count(i) + 1
               forest%oob_prediction(i) = forest%oob_prediction(i) + &
                  (tree_pred(i) - forest%oob_prediction(i)) / real(forest%oob_count(i), dp)
            end if
         end do
         n_oob = count(forest%oob_count > 0)
         if (n_oob > 0) then
            if (opt%bias_correct) then
               call simple_linear_bias(forest%oob_prediction, y, forest%oob_count > 0, intercept, slope)
               mse = sum((y - (intercept + slope * forest%oob_prediction)) ** 2, mask=forest%oob_count > 0) / real(n_oob, dp)
               forest%bias_intercept = intercept
               forest%bias_slope = slope
            else
               mse = sum((y - forest%oob_prediction) ** 2, mask=forest%oob_count > 0) / real(n_oob, dp)
            end if
            forest%mse_curve(t) = mse
         end if

         if (opt%proximity) call accumulate_proximity(leaf, inbag_count, opt%oob_proximity, prox_count, prox_denom)
         if (opt%importance) then
            call tree_var_used(forest%trees(t), used)
            call regression_tree_importance(forest%trees(t), x, y, cats, inbag_count, tree_pred, used, &
               max(1, opt%n_perm), rng, imp_sum, imp_sumsq)
         end if
      end do

      forest%importance_gini = forest%importance_gini / real(ntree, dp)
      if (opt%importance) call finish_importance(imp_sum, imp_sumsq, ntree, forest%importance_accuracy, forest%importance_sd)
      if (opt%proximity) call finish_proximity(prox_count, prox_denom, ntree, opt%oob_proximity, forest%proximity)
   end subroutine fit_regression

   subroutine predict_regression(forest, x, prediction, per_tree, terminal_nodes)
      type(rf_regression_forest), intent(in) :: forest
      real(dp), intent(in) :: x(:,:)
      real(dp), intent(out) :: prediction(:)
      real(dp), intent(out), optional :: per_tree(:,:)
      integer, intent(out), optional :: terminal_nodes(:,:)
      real(dp), allocatable :: tree_pred(:)
      integer, allocatable :: leaf(:)
      integer :: n, t

      n = size(x, 1)
      if (size(prediction) /= n) error stop 'predict_regression: prediction has wrong size'
      if (size(x, 2) /= forest%nvar) error stop 'predict_regression: x has wrong number of variables'
      if (present(per_tree)) then
         if (size(per_tree, 1) /= n .or. size(per_tree, 2) /= size(forest%trees)) &
            error stop 'predict_regression: per_tree has wrong shape'
      end if
      if (present(terminal_nodes)) then
         if (size(terminal_nodes, 1) /= n .or. size(terminal_nodes, 2) /= size(forest%trees)) &
            error stop 'predict_regression: terminal_nodes has wrong shape'
      end if

      allocate(tree_pred(n), leaf(n))
      prediction = 0.0_dp
      do t = 1, size(forest%trees)
         call predict_reg_tree(forest%trees(t), x, forest%ncat, tree_pred, leaf)
         prediction = prediction + tree_pred
         if (present(per_tree)) per_tree(:, t) = tree_pred
         if (present(terminal_nodes)) terminal_nodes(:, t) = leaf
      end do
      prediction = prediction / real(size(forest%trees), dp)
      prediction = forest%bias_intercept + forest%bias_slope * prediction
   end subroutine predict_regression

   subroutine regression_tree_importance(tree, x, y, ncat, inbag, base_pred, used, nperm, rng, imp_sum, imp_sumsq)
      type(rf_tree), intent(in) :: tree
      real(dp), intent(in) :: x(:,:), y(:), base_pred(:)
      integer, intent(in) :: ncat(:), inbag(:), nperm
      logical, intent(in) :: used(:)
      type(rf_rng_state), intent(inout) :: rng
      real(dp), intent(inout) :: imp_sum(:), imp_sumsq(:)
      real(dp), allocatable :: xperm(:,:), vals(:), pred(:)
      integer, allocatable :: oob_idx(:)
      integer :: m, k, i, noob
      real(dp) :: base_mse, perm_mse, delta

      noob = count(inbag == 0)
      if (noob <= 0) return
      allocate(oob_idx(noob), xperm(size(x, 1), size(x, 2)), vals(noob), pred(size(y)))
      oob_idx = pack([(i, i = 1, size(y))], inbag == 0)
      base_mse = sum((base_pred(oob_idx) - y(oob_idx)) ** 2) / real(noob, dp)

      do m = 1, size(used)
         if (.not. used(m)) cycle
         perm_mse = 0.0_dp
         do k = 1, nperm
            xperm = x
            vals = x(oob_idx, m)
            call shuffle_real(rng, vals)
            xperm(oob_idx, m) = vals
            call predict_reg_tree(tree, xperm, ncat, pred)
            perm_mse = perm_mse + sum((pred(oob_idx) - y(oob_idx)) ** 2) / real(noob, dp)
         end do
         perm_mse = perm_mse / real(nperm, dp)
         delta = perm_mse - base_mse
         imp_sum(m) = imp_sum(m) + delta
         imp_sumsq(m) = imp_sumsq(m) + delta * delta
      end do
   end subroutine regression_tree_importance

   subroutine simple_linear_bias(x, y, mask, intercept, slope)
      real(dp), intent(in) :: x(:), y(:)
      logical, intent(in) :: mask(:)
      real(dp), intent(out) :: intercept, slope
      real(dp) :: mx, my, sxx, sxy
      integer :: n

      n = count(mask)
      if (n <= 1) then
         intercept = 0.0_dp
         slope = 1.0_dp
         return
      end if
      mx = sum(x, mask=mask) / real(n, dp)
      my = sum(y, mask=mask) / real(n, dp)
      sxx = sum((x - mx) ** 2, mask=mask)
      sxy = sum((x - mx) * (y - my), mask=mask)
      if (sxx <= epsilon(1.0_dp)) then
         slope = 0.0_dp
         intercept = my
      else
         slope = sxy / sxx
         intercept = my - slope * mx
      end if
   end subroutine simple_linear_bias

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

   subroutine finish_importance(sum_imp, sumsq_imp, ntree, mean_imp, sd_imp)
      real(dp), intent(in) :: sum_imp(:), sumsq_imp(:)
      integer, intent(in) :: ntree
      real(dp), intent(out) :: mean_imp(:), sd_imp(:)
      real(dp) :: v
      integer :: i

      mean_imp = sum_imp / real(ntree, dp)
      do i = 1, size(mean_imp)
         v = max(0.0_dp, sumsq_imp(i) / real(ntree, dp) - mean_imp(i) ** 2)
         sd_imp(i) = sqrt(v / real(ntree, dp))
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

end module rf_regression
