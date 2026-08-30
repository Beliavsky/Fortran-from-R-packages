! SPDX-License-Identifier: GPL-2.0-or-later
! Modern Fortran translation of gbm3 3.0.3; see NOTICE.md and upstream/.
module gbm3_core
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
   use gbm3_kinds, only : dp
   use gbm3_constants
   use gbm3_math, only : quiet_nan, sigmoid
   use gbm3_types, only : gbm_options, gbm_model, gbm_tree
   use gbm3_tree, only : grow_tree, adjust_tree, predict_tree, tree_relative_influence
   use gbm3_distributions, only : dist_init_f, dist_working_response, dist_deviance, &
                                  dist_fit_best_constant, dist_bag_improvement
   use gbm3_pairwise, only : pairwise_make_bag, pairwise_working_response, pairwise_deviance, &
                             pairwise_fit_best_constant, pairwise_bag_improvement
   use gbm3_cox, only : cox_working_response, cox_deviance, cox_fit_best_constant, cox_bag_improvement
   implicit none
   private

   interface gbm_fit
      module procedure gbm_fit_vector
      module procedure gbm_fit_cox
   end interface gbm_fit

   interface gbm_continue
      module procedure gbm_continue_vector
      module procedure gbm_continue_cox
   end interface gbm_continue

   public :: gbm_fit, gbm_fit_vector, gbm_fit_cox
   public :: gbm_continue, gbm_continue_vector, gbm_continue_cox
   public :: gbm_predict, gbm_predict_response, gbm_predict_trees, gbm_predict_staged, gbm_predict_response_staged
   public :: gbm_partial_dependence, gbm_relative_influence
   public :: gbm_set_seed

contains

   subroutine gbm_fit_vector(x, y, model, options, weight, offset, var_classes, monotone, id, group)
      real(dp), intent(in) :: x(:, :), y(:)
      type(gbm_model), intent(out) :: model
      type(gbm_options), intent(in), optional :: options
      real(dp), intent(in), optional :: weight(:), offset(:)
      integer, intent(in), optional :: var_classes(:), monotone(:), id(:), group(:)

      type(gbm_options) :: opt
      real(dp), allocatable :: w(:), off(:), f(:), residual(:), hessian(:), delta(:)
      integer, allocatable :: vc(:), mono(:), ids(:), grp(:), assignment(:)
      logical, allocatable :: in_bag(:)
      integer :: n, p, ntrain, ntree, i

      n = size(x, 1)
      p = size(x, 2)
      if (size(y) /= n) error stop "gbm_fit: response length does not match x"
      if (n <= 0 .or. p <= 0) error stop "gbm_fit: x must be non-empty"

      opt = gbm_options()
      if (present(options)) opt = options
      call validate_options(opt, p)
      if (opt%distribution == GBM_COXPH) error stop "gbm_fit: Cox PH requires a survival response matrix"
      if (opt%distribution == GBM_PAIRWISE .and. .not. present(group)) &
         error stop "gbm_fit: pairwise distribution requires group"

      ntrain = opt%num_train
      if (ntrain <= 0) ntrain = n
      if (ntrain > n) error stop "gbm_fit: num_train exceeds number of rows"
      if (ntrain < 1) error stop "gbm_fit: num_train must be positive"

      allocate(w(n), off(n), f(n), residual(ntrain), delta(ntrain), in_bag(ntrain))
      allocate(vc(p), mono(p), ids(ntrain))
      w = 1.0_dp
      if (present(weight)) then
         if (size(weight) /= n) error stop "gbm_fit: weight length mismatch"
         w = weight
      end if
      if (any(w < 0.0_dp)) error stop "gbm_fit: weights must be nonnegative"
      if (sum(w(1:ntrain)) <= 0.0_dp) error stop "gbm_fit: positive training weight is required"

      off = 0.0_dp
      if (present(offset)) then
         if (size(offset) /= n) error stop "gbm_fit: offset length mismatch"
         off = offset
      end if

      vc = 0
      if (present(var_classes)) then
         if (size(var_classes) /= p) error stop "gbm_fit: var_classes length mismatch"
         vc = var_classes
      end if
      if (any(vc < 0)) error stop "gbm_fit: var_classes values must be nonnegative"

      mono = 0
      if (present(monotone)) then
         if (size(monotone) /= p) error stop "gbm_fit: monotone length mismatch"
         mono = monotone
      end if
      if (any(abs(mono) > 1)) error stop "gbm_fit: monotone values must be -1, 0, or 1"
      where (vc > 0) mono = 0

      do i = 1, ntrain
         ids(i) = i
      end do
      if (present(id)) then
         if (size(id) /= n) error stop "gbm_fit: id length mismatch"
         ids = id(1:ntrain)
      end if

      if (opt%distribution == GBM_PAIRWISE) then
         allocate(grp(n))
         if (size(group) /= n) error stop "gbm_fit: group length mismatch"
         grp = group
         call validate_groups(grp)
         if (ntrain < n) then
            if (grp(ntrain) == grp(ntrain + 1)) error stop "gbm_fit: num_train cannot split a pairwise group"
         end if
         model%init_f = 0.0_dp
      else
         model%init_f = dist_init_f(y(1:ntrain), off(1:ntrain), w(1:ntrain), opt)
      end if
      f = model%init_f

      call initialize_model(model, opt, p, vc, mono)
      model%n_rows = n
      model%n_train = ntrain
      model%init_f = merge(0.0_dp, model%init_f, opt%distribution == GBM_PAIRWISE)
      if (opt%distribution == GBM_PAIRWISE) model%init_f = 0.0_dp
      f = model%init_f

      if (opt%distribution == GBM_PAIRWISE) allocate(hessian(ntrain))

      do ntree = 1, opt%num_trees
         if (opt%distribution == GBM_PAIRWISE) then
            call pairwise_make_bag(grp(1:ntrain), opt%bag_fraction, in_bag)
            call pairwise_working_response(y(1:ntrain), grp(1:ntrain), off(1:ntrain), f(1:ntrain), &
                                           w(1:ntrain), in_bag, opt, residual, hessian)
         else
            call make_bag(ids, opt%bag_fraction, in_bag)
            call dist_working_response(y(1:ntrain), off(1:ntrain), f(1:ntrain), opt, residual)
         end if

         call grow_tree(x(1:ntrain, :), residual, w(1:ntrain), in_bag, vc, mono, &
                        opt%interaction_depth, opt%min_num_obs_in_node, opt%num_features, &
                        opt%shrinkage, model%trees(ntree), assignment)

         if (opt%distribution == GBM_PAIRWISE) then
            call pairwise_fit_best_constant(w(1:ntrain), residual, hessian, in_bag, assignment, model%trees(ntree))
         else
            call dist_fit_best_constant(y(1:ntrain), off(1:ntrain), w(1:ntrain), f(1:ntrain), residual, &
                                        in_bag, assignment, opt%min_num_obs_in_node, opt, model%trees(ntree))
         end if

         call adjust_tree(model%trees(ntree), opt%min_num_obs_in_node, assignment, delta)

         if (opt%distribution == GBM_PAIRWISE) then
            model%oob_improvement(ntree) = pairwise_bag_improvement(y(1:ntrain), grp(1:ntrain), &
                 off(1:ntrain), f(1:ntrain), w(1:ntrain), delta, in_bag, opt)
         else
            model%oob_improvement(ntree) = dist_bag_improvement(y(1:ntrain), off(1:ntrain), &
                 w(1:ntrain), f(1:ntrain), delta, in_bag, opt)
         end if

         f(1:ntrain) = f(1:ntrain) + opt%shrinkage * delta
         if (opt%distribution == GBM_PAIRWISE) then
            model%train_error(ntree) = pairwise_deviance(y(1:ntrain), grp(1:ntrain), off(1:ntrain), &
                                                         f(1:ntrain), w(1:ntrain), opt)
         else
            model%train_error(ntree) = dist_deviance(y(1:ntrain), off(1:ntrain), w(1:ntrain), f(1:ntrain), opt)
         end if

         if (ntrain < n) then
            do i = ntrain + 1, n
               f(i) = f(i) + predict_tree(model%trees(ntree), x(i, :))
            end do
            if (opt%distribution == GBM_PAIRWISE) then
               model%validation_error(ntree) = pairwise_deviance(y(ntrain + 1:n), grp(ntrain + 1:n), &
                    off(ntrain + 1:n), f(ntrain + 1:n), w(ntrain + 1:n), opt)
            else
               model%validation_error(ntree) = dist_deviance(y(ntrain + 1:n), off(ntrain + 1:n), &
                    w(ntrain + 1:n), f(ntrain + 1:n), opt)
            end if
         else
            model%validation_error(ntree) = quiet_nan()
         end if
      end do
      model%n_trees = opt%num_trees
      model%fitted = f
   end subroutine gbm_fit_vector

   subroutine gbm_fit_cox(x, surv, model, options, weight, offset, var_classes, monotone, id, strata)
      real(dp), intent(in) :: x(:, :), surv(:, :)
      type(gbm_model), intent(out) :: model
      type(gbm_options), intent(in), optional :: options
      real(dp), intent(in), optional :: weight(:), offset(:)
      integer, intent(in), optional :: var_classes(:), monotone(:), id(:), strata(:)

      type(gbm_options) :: opt
      real(dp), allocatable :: w(:), off(:), f(:), residual(:), delta(:)
      integer, allocatable :: vc(:), mono(:), ids(:), str(:), assignment(:)
      logical, allocatable :: in_bag(:)
      integer :: n, p, ntrain, ntree, i

      n = size(x, 1)
      p = size(x, 2)
      if (size(surv, 1) /= n) error stop "gbm_fit: survival response row mismatch"
      if (size(surv, 2) /= 2 .and. size(surv, 2) /= 3) &
         error stop "gbm_fit: survival response must have (time,status) or (start,stop,status) columns"
      if (n <= 0 .or. p <= 0) error stop "gbm_fit: x must be non-empty"

      opt = gbm_options(distribution=GBM_COXPH)
      if (present(options)) opt = options
      if (opt%distribution /= GBM_COXPH) error stop "gbm_fit: survival response requires GBM_COXPH"
      call validate_options(opt, p)

      ntrain = opt%num_train
      if (ntrain <= 0) ntrain = n
      if (ntrain > n) error stop "gbm_fit: num_train exceeds number of rows"

      allocate(w(n), off(n), f(n), residual(ntrain), delta(ntrain), in_bag(ntrain))
      allocate(vc(p), mono(p), ids(ntrain), str(n), assignment(ntrain))
      w = 1.0_dp
      if (present(weight)) then
         if (size(weight) /= n) error stop "gbm_fit: weight length mismatch"
         w = weight
      end if
      if (any(w < 0.0_dp)) error stop "gbm_fit: weights must be nonnegative"
      if (sum(w(1:ntrain)) <= 0.0_dp) error stop "gbm_fit: positive training weight is required"

      off = 0.0_dp
      if (present(offset)) then
         if (size(offset) /= n) error stop "gbm_fit: offset length mismatch"
         off = offset
      end if

      vc = 0
      if (present(var_classes)) then
         if (size(var_classes) /= p) error stop "gbm_fit: var_classes length mismatch"
         vc = var_classes
      end if
      if (any(vc < 0)) error stop "gbm_fit: var_classes values must be nonnegative"

      mono = 0
      if (present(monotone)) then
         if (size(monotone) /= p) error stop "gbm_fit: monotone length mismatch"
         mono = monotone
      end if
      if (any(abs(mono) > 1)) error stop "gbm_fit: monotone values must be -1, 0, or 1"
      where (vc > 0) mono = 0

      do i = 1, ntrain
         ids(i) = i
      end do
      if (present(id)) then
         if (size(id) /= n) error stop "gbm_fit: id length mismatch"
         ids = id(1:ntrain)
      end if

      str = 1
      if (present(strata)) then
         if (size(strata) /= n) error stop "gbm_fit: strata length mismatch"
         str = strata
      end if

      model%init_f = 0.0_dp
      f = 0.0_dp
      call initialize_model(model, opt, p, vc, mono)
      model%n_rows = n
      model%n_train = ntrain
      model%init_f = 0.0_dp

      do ntree = 1, opt%num_trees
         call make_bag(ids, opt%bag_fraction, in_bag)
         call cox_working_response(surv(1:ntrain, :), str(1:ntrain), off(1:ntrain), f(1:ntrain), &
                                   w(1:ntrain), in_bag, opt, residual)

         call grow_tree(x(1:ntrain, :), residual, w(1:ntrain), in_bag, vc, mono, &
                        opt%interaction_depth, opt%min_num_obs_in_node, opt%num_features, &
                        opt%shrinkage, model%trees(ntree), assignment)
         call cox_fit_best_constant(surv(1:ntrain, :), str(1:ntrain), off(1:ntrain), f(1:ntrain), &
                                    w(1:ntrain), in_bag, assignment, opt%min_num_obs_in_node, opt, &
                                    model%trees(ntree))
         call adjust_tree(model%trees(ntree), opt%min_num_obs_in_node, assignment, delta)
         model%oob_improvement(ntree) = cox_bag_improvement(surv(1:ntrain, :), str(1:ntrain), &
              off(1:ntrain), f(1:ntrain), w(1:ntrain), delta, in_bag, opt)

         f(1:ntrain) = f(1:ntrain) + opt%shrinkage * delta
         model%train_error(ntree) = cox_deviance(surv(1:ntrain, :), str(1:ntrain), off(1:ntrain), &
                                                 f(1:ntrain), w(1:ntrain), opt)

         if (ntrain < n) then
            do i = ntrain + 1, n
               f(i) = f(i) + predict_tree(model%trees(ntree), x(i, :))
            end do
            model%validation_error(ntree) = cox_deviance(surv(ntrain + 1:n, :), str(ntrain + 1:n), &
                 off(ntrain + 1:n), f(ntrain + 1:n), w(ntrain + 1:n), opt)
         else
            model%validation_error(ntree) = quiet_nan()
         end if
      end do
      model%n_trees = opt%num_trees
      model%fitted = f
   end subroutine gbm_fit_cox

   subroutine gbm_continue_vector(model, x, y, additional_trees, weight, offset, id, group)
      type(gbm_model), intent(inout) :: model
      real(dp), intent(in) :: x(:, :), y(:)
      integer, intent(in) :: additional_trees
      real(dp), intent(in), optional :: weight(:), offset(:)
      integer, intent(in), optional :: id(:), group(:)

      type(gbm_options) :: opt
      real(dp), allocatable :: w(:), off(:), f(:), residual(:), hessian(:), delta(:)
      integer, allocatable :: ids(:), grp(:), assignment(:)
      logical, allocatable :: in_bag(:)
      integer :: n, p, ntrain, old_ntrees, total_trees, ntree, i

      if (additional_trees < 1) error stop "gbm_continue: additional_trees must be positive"
      n = size(x, 1)
      p = size(x, 2)
      if (size(y) /= n) error stop "gbm_continue: response length does not match x"
      if (p /= model%n_features) error stop "gbm_continue: feature count differs from original fit"
      if (model%n_rows > 0 .and. n /= model%n_rows) &
         error stop "gbm_continue: row count differs from original fit"
      if (model%n_trees < 1 .or. .not. allocated(model%trees)) &
         error stop "gbm_continue: model has no fitted trees"

      opt = model%options
      if (opt%distribution == GBM_COXPH) error stop "gbm_continue: Cox PH requires a survival response matrix"
      if (opt%distribution == GBM_PAIRWISE .and. .not. present(group)) &
         error stop "gbm_continue: pairwise distribution requires group"
      old_ntrees = model%n_trees
      total_trees = old_ntrees + additional_trees
      opt%num_trees = total_trees
      call validate_options(opt, p)

      ntrain = model%n_train
      if (ntrain <= 0) then
         ntrain = model%options%num_train
         if (ntrain <= 0) ntrain = n
      end if
      if (ntrain > n .or. ntrain < 1) error stop "gbm_continue: invalid stored training row count"

      allocate(w(n), off(n), f(n), residual(ntrain), delta(ntrain), in_bag(ntrain), ids(ntrain))
      w = 1.0_dp
      if (present(weight)) then
         if (size(weight) /= n) error stop "gbm_continue: weight length mismatch"
         w = weight
      end if
      if (any(w < 0.0_dp)) error stop "gbm_continue: weights must be nonnegative"
      if (sum(w(1:ntrain)) <= 0.0_dp) error stop "gbm_continue: positive training weight is required"

      off = 0.0_dp
      if (present(offset)) then
         if (size(offset) /= n) error stop "gbm_continue: offset length mismatch"
         off = offset
      end if

      do i = 1, ntrain
         ids(i) = i
      end do
      if (present(id)) then
         if (size(id) /= n) error stop "gbm_continue: id length mismatch"
         ids = id(1:ntrain)
      end if

      if (opt%distribution == GBM_PAIRWISE) then
         allocate(grp(n), hessian(ntrain))
         if (size(group) /= n) error stop "gbm_continue: group length mismatch"
         grp = group
         call validate_groups(grp)
         if (ntrain < n) then
            if (grp(ntrain) == grp(ntrain + 1)) error stop "gbm_continue: stored num_train splits a pairwise group"
         end if
      end if

      if (allocated(model%fitted) .and. size(model%fitted) == n) then
         f = model%fitted
      else
         f = gbm_predict(model, x)
      end if
      call append_model_storage(model, total_trees)

      do ntree = old_ntrees + 1, total_trees
         if (opt%distribution == GBM_PAIRWISE) then
            call pairwise_make_bag(grp(1:ntrain), opt%bag_fraction, in_bag)
            call pairwise_working_response(y(1:ntrain), grp(1:ntrain), off(1:ntrain), f(1:ntrain), &
                                           w(1:ntrain), in_bag, opt, residual, hessian)
         else
            call make_bag(ids, opt%bag_fraction, in_bag)
            call dist_working_response(y(1:ntrain), off(1:ntrain), f(1:ntrain), opt, residual)
         end if

         call grow_tree(x(1:ntrain, :), residual, w(1:ntrain), in_bag, model%var_classes, model%monotone, &
                        opt%interaction_depth, opt%min_num_obs_in_node, opt%num_features, &
                        opt%shrinkage, model%trees(ntree), assignment)

         if (opt%distribution == GBM_PAIRWISE) then
            call pairwise_fit_best_constant(w(1:ntrain), residual, hessian, in_bag, assignment, model%trees(ntree))
         else
            call dist_fit_best_constant(y(1:ntrain), off(1:ntrain), w(1:ntrain), f(1:ntrain), residual, &
                                        in_bag, assignment, opt%min_num_obs_in_node, opt, model%trees(ntree))
         end if

         call adjust_tree(model%trees(ntree), opt%min_num_obs_in_node, assignment, delta)
         if (opt%distribution == GBM_PAIRWISE) then
            model%oob_improvement(ntree) = pairwise_bag_improvement(y(1:ntrain), grp(1:ntrain), &
                 off(1:ntrain), f(1:ntrain), w(1:ntrain), delta, in_bag, opt)
         else
            model%oob_improvement(ntree) = dist_bag_improvement(y(1:ntrain), off(1:ntrain), &
                 w(1:ntrain), f(1:ntrain), delta, in_bag, opt)
         end if

         f(1:ntrain) = f(1:ntrain) + opt%shrinkage * delta
         if (opt%distribution == GBM_PAIRWISE) then
            model%train_error(ntree) = pairwise_deviance(y(1:ntrain), grp(1:ntrain), off(1:ntrain), &
                                                         f(1:ntrain), w(1:ntrain), opt)
         else
            model%train_error(ntree) = dist_deviance(y(1:ntrain), off(1:ntrain), w(1:ntrain), f(1:ntrain), opt)
         end if

         if (ntrain < n) then
            do i = ntrain + 1, n
               f(i) = f(i) + predict_tree(model%trees(ntree), x(i, :))
            end do
            if (opt%distribution == GBM_PAIRWISE) then
               model%validation_error(ntree) = pairwise_deviance(y(ntrain + 1:n), grp(ntrain + 1:n), &
                    off(ntrain + 1:n), f(ntrain + 1:n), w(ntrain + 1:n), opt)
            else
               model%validation_error(ntree) = dist_deviance(y(ntrain + 1:n), off(ntrain + 1:n), &
                    w(ntrain + 1:n), f(ntrain + 1:n), opt)
            end if
         else
            model%validation_error(ntree) = quiet_nan()
         end if
      end do

      model%n_trees = total_trees
      model%options%num_trees = total_trees
      model%fitted = f
   end subroutine gbm_continue_vector

   subroutine gbm_continue_cox(model, x, surv, additional_trees, weight, offset, id, strata)
      type(gbm_model), intent(inout) :: model
      real(dp), intent(in) :: x(:, :), surv(:, :)
      integer, intent(in) :: additional_trees
      real(dp), intent(in), optional :: weight(:), offset(:)
      integer, intent(in), optional :: id(:), strata(:)

      type(gbm_options) :: opt
      real(dp), allocatable :: w(:), off(:), f(:), residual(:), delta(:)
      integer, allocatable :: ids(:), str(:), assignment(:)
      logical, allocatable :: in_bag(:)
      integer :: n, p, ntrain, old_ntrees, total_trees, ntree, i

      if (additional_trees < 1) error stop "gbm_continue: additional_trees must be positive"
      n = size(x, 1)
      p = size(x, 2)
      if (size(surv, 1) /= n) error stop "gbm_continue: survival response row mismatch"
      if (size(surv, 2) /= 2 .and. size(surv, 2) /= 3) &
         error stop "gbm_continue: survival response must have (time,status) or (start,stop,status) columns"
      if (p /= model%n_features) error stop "gbm_continue: feature count differs from original fit"
      if (model%n_rows > 0 .and. n /= model%n_rows) &
         error stop "gbm_continue: row count differs from original fit"
      if (model%n_trees < 1 .or. .not. allocated(model%trees)) &
         error stop "gbm_continue: model has no fitted trees"

      opt = model%options
      if (opt%distribution /= GBM_COXPH) error stop "gbm_continue: model is not Cox PH"
      old_ntrees = model%n_trees
      total_trees = old_ntrees + additional_trees
      opt%num_trees = total_trees
      call validate_options(opt, p)

      ntrain = model%n_train
      if (ntrain <= 0) then
         ntrain = model%options%num_train
         if (ntrain <= 0) ntrain = n
      end if
      if (ntrain > n .or. ntrain < 1) error stop "gbm_continue: invalid stored training row count"

      allocate(w(n), off(n), f(n), residual(ntrain), delta(ntrain), in_bag(ntrain), ids(ntrain), str(n))
      w = 1.0_dp
      if (present(weight)) then
         if (size(weight) /= n) error stop "gbm_continue: weight length mismatch"
         w = weight
      end if
      if (any(w < 0.0_dp)) error stop "gbm_continue: weights must be nonnegative"
      if (sum(w(1:ntrain)) <= 0.0_dp) error stop "gbm_continue: positive training weight is required"

      off = 0.0_dp
      if (present(offset)) then
         if (size(offset) /= n) error stop "gbm_continue: offset length mismatch"
         off = offset
      end if
      do i = 1, ntrain
         ids(i) = i
      end do
      if (present(id)) then
         if (size(id) /= n) error stop "gbm_continue: id length mismatch"
         ids = id(1:ntrain)
      end if
      str = 1
      if (present(strata)) then
         if (size(strata) /= n) error stop "gbm_continue: strata length mismatch"
         str = strata
      end if

      if (allocated(model%fitted) .and. size(model%fitted) == n) then
         f = model%fitted
      else
         f = gbm_predict(model, x)
      end if
      call append_model_storage(model, total_trees)

      do ntree = old_ntrees + 1, total_trees
         call make_bag(ids, opt%bag_fraction, in_bag)
         call cox_working_response(surv(1:ntrain, :), str(1:ntrain), off(1:ntrain), f(1:ntrain), &
                                   w(1:ntrain), in_bag, opt, residual)
         call grow_tree(x(1:ntrain, :), residual, w(1:ntrain), in_bag, model%var_classes, model%monotone, &
                        opt%interaction_depth, opt%min_num_obs_in_node, opt%num_features, &
                        opt%shrinkage, model%trees(ntree), assignment)
         call cox_fit_best_constant(surv(1:ntrain, :), str(1:ntrain), off(1:ntrain), f(1:ntrain), &
                                    w(1:ntrain), in_bag, assignment, opt%min_num_obs_in_node, opt, &
                                    model%trees(ntree))
         call adjust_tree(model%trees(ntree), opt%min_num_obs_in_node, assignment, delta)
         model%oob_improvement(ntree) = cox_bag_improvement(surv(1:ntrain, :), str(1:ntrain), &
              off(1:ntrain), f(1:ntrain), w(1:ntrain), delta, in_bag, opt)

         f(1:ntrain) = f(1:ntrain) + opt%shrinkage * delta
         model%train_error(ntree) = cox_deviance(surv(1:ntrain, :), str(1:ntrain), off(1:ntrain), &
                                                 f(1:ntrain), w(1:ntrain), opt)
         if (ntrain < n) then
            do i = ntrain + 1, n
               f(i) = f(i) + predict_tree(model%trees(ntree), x(i, :))
            end do
            model%validation_error(ntree) = cox_deviance(surv(ntrain + 1:n, :), str(ntrain + 1:n), &
                 off(ntrain + 1:n), f(ntrain + 1:n), w(ntrain + 1:n), opt)
         else
            model%validation_error(ntree) = quiet_nan()
         end if
      end do

      model%n_trees = total_trees
      model%options%num_trees = total_trees
      model%fitted = f
   end subroutine gbm_continue_cox

   subroutine append_model_storage(model, new_total)
      type(gbm_model), intent(inout) :: model
      integer, intent(in) :: new_total
      type(gbm_tree), allocatable :: old_trees(:)
      real(dp), allocatable :: old_train(:), old_validation(:), old_oob(:)
      integer :: old_n

      old_n = model%n_trees
      if (new_total <= old_n) error stop "append_model_storage: new size must exceed current tree count"
      if (.not. allocated(model%trees) .or. .not. allocated(model%train_error) .or. &
          .not. allocated(model%validation_error) .or. .not. allocated(model%oob_improvement)) &
         error stop "append_model_storage: incomplete model storage"

      allocate(old_trees(old_n), old_train(old_n), old_validation(old_n), old_oob(old_n))
      if (old_n > 0) then
         old_trees = model%trees(1:old_n)
         old_train = model%train_error(1:old_n)
         old_validation = model%validation_error(1:old_n)
         old_oob = model%oob_improvement(1:old_n)
      end if
      deallocate(model%trees, model%train_error, model%validation_error, model%oob_improvement)
      allocate(model%trees(new_total), model%train_error(new_total), &
               model%validation_error(new_total), model%oob_improvement(new_total))
      model%train_error = quiet_nan()
      model%validation_error = quiet_nan()
      model%oob_improvement = quiet_nan()
      if (old_n > 0) then
         model%trees(1:old_n) = old_trees(1:old_n)
         model%train_error(1:old_n) = old_train(1:old_n)
         model%validation_error(1:old_n) = old_validation(1:old_n)
         model%oob_improvement(1:old_n) = old_oob(1:old_n)
      end if
   end subroutine append_model_storage

   subroutine initialize_model(model, options, p, var_classes, monotone)
      type(gbm_model), intent(inout) :: model
      type(gbm_options), intent(in) :: options
      integer, intent(in) :: p, var_classes(:), monotone(:)
      model%options = options
      model%n_features = p
      model%n_rows = 0
      model%n_train = 0
      model%n_trees = 0
      if (allocated(model%trees)) deallocate(model%trees)
      if (allocated(model%train_error)) deallocate(model%train_error)
      if (allocated(model%validation_error)) deallocate(model%validation_error)
      if (allocated(model%oob_improvement)) deallocate(model%oob_improvement)
      if (allocated(model%fitted)) deallocate(model%fitted)
      if (allocated(model%var_classes)) deallocate(model%var_classes)
      if (allocated(model%monotone)) deallocate(model%monotone)
      allocate(model%trees(options%num_trees), model%train_error(options%num_trees), &
               model%validation_error(options%num_trees), model%oob_improvement(options%num_trees), &
               model%var_classes(p), model%monotone(p))
      model%train_error = quiet_nan()
      model%validation_error = quiet_nan()
      model%oob_improvement = quiet_nan()
      model%var_classes = var_classes
      model%monotone = monotone
   end subroutine initialize_model

   subroutine validate_options(options, p)
      type(gbm_options), intent(in) :: options
      integer, intent(in) :: p
      if (options%distribution < GBM_GAUSSIAN .or. options%distribution > GBM_PAIRWISE) &
         error stop "gbm_fit: unknown distribution"
      if (options%num_trees < 1) error stop "gbm_fit: num_trees must be positive"
      if (options%interaction_depth < 1) error stop "gbm_fit: interaction_depth must be positive"
      if (options%min_num_obs_in_node < 1) error stop "gbm_fit: min_num_obs_in_node must be positive"
      if (options%shrinkage <= 0.0_dp) error stop "gbm_fit: shrinkage must be positive"
      if (options%bag_fraction <= 0.0_dp .or. options%bag_fraction > 1.0_dp) &
         error stop "gbm_fit: bag_fraction must be in (0,1]"
      if (options%num_features < 0 .or. options%num_features > p) &
         error stop "gbm_fit: num_features must be 0..number of predictors"
      if (options%quantile_alpha <= 0.0_dp .or. options%quantile_alpha >= 1.0_dp) &
         error stop "gbm_fit: quantile_alpha must be in (0,1)"
      if (options%t_df <= 0.0_dp) error stop "gbm_fit: t_df must be positive"
      if (options%tweedie_power <= 1.0_dp .or. options%tweedie_power >= 2.0_dp) &
         error stop "gbm_fit: tweedie_power must be in (1,2)"
      if (options%cox_prior_node_coeff_var <= 0.0_dp) &
         error stop "gbm_fit: cox_prior_node_coeff_var must be positive"
   end subroutine validate_options

   subroutine validate_groups(group)
      integer, intent(in) :: group(:)
      integer :: i, j
      ! gbm3's pairwise implementation requires observations belonging to a
      ! group to be contiguous. Reject a group that reappears later.
      i = 1
      do while (i <= size(group))
         j = i + 1
         do while (j <= size(group))
            if (group(j) /= group(i)) exit
            j = j + 1
         end do
         if (j <= size(group)) then
            if (any(group(j:) == group(i))) error stop "gbm_fit: pairwise groups must be contiguous"
         end if
         i = j
      end do
   end subroutine validate_groups

   subroutine make_bag(id, bag_fraction, in_bag)
      integer, intent(in) :: id(:)
      real(dp), intent(in) :: bag_fraction
      logical, intent(out) :: in_bag(size(id))
      integer, allocatable :: unique_id(:)
      integer :: i, j, n_unique, target, seen, chosen
      logical :: exists
      real(dp) :: u

      in_bag = .false.
      if (size(id) == 0) return
      allocate(unique_id(size(id)))
      n_unique = 0
      do i = 1, size(id)
         exists = .false.
         do j = 1, n_unique
            if (id(i) == unique_id(j)) then
               exists = .true.
               exit
            end if
         end do
         if (.not. exists) then
            n_unique = n_unique + 1
            unique_id(n_unique) = id(i)
         end if
      end do

      target = int(bag_fraction * real(n_unique, dp))
      if (target <= 0) target = 1
      seen = 0
      chosen = 0
      do i = 1, n_unique
         if (chosen >= target) exit
         call random_number(u)
         if (u * real(n_unique - seen, dp) < real(target - chosen, dp)) then
            where (id == unique_id(i)) in_bag = .true.
            chosen = chosen + 1
         end if
         seen = seen + 1
      end do
   end subroutine make_bag

   function gbm_predict(model, x, offset, n_trees) result(pred)
      type(gbm_model), intent(in) :: model
      real(dp), intent(in) :: x(:, :)
      real(dp), intent(in), optional :: offset(:)
      integer, intent(in), optional :: n_trees
      real(dp), allocatable :: pred(:)
      integer :: n, nt, i, t

      n = size(x, 1)
      if (size(x, 2) /= model%n_features) error stop "gbm_predict: feature count mismatch"
      nt = model%n_trees
      if (present(n_trees)) nt = n_trees
      if (nt < 0 .or. nt > model%n_trees) error stop "gbm_predict: invalid n_trees"

      allocate(pred(n))
      pred = model%init_f
      if (present(offset)) then
         if (size(offset) /= n) error stop "gbm_predict: offset length mismatch"
         pred = pred + offset
      end if
      do t = 1, nt
         do i = 1, n
            pred(i) = pred(i) + predict_tree(model%trees(t), x(i, :))
         end do
      end do
   end function gbm_predict

   function gbm_predict_response(model, x, offset, n_trees) result(pred)
      type(gbm_model), intent(in) :: model
      real(dp), intent(in) :: x(:, :)
      real(dp), intent(in), optional :: offset(:)
      integer, intent(in), optional :: n_trees
      real(dp), allocatable :: pred(:)
      integer :: i

      pred = gbm_predict(model, x, offset, n_trees)
      select case (model%options%distribution)
      case (GBM_BERNOULLI)
         do i = 1, size(pred)
            pred(i) = sigmoid(pred(i))
         end do
      case (GBM_ADABOOST)
         do i = 1, size(pred)
            pred(i) = sigmoid(2.0_dp * pred(i))
         end do
      case (GBM_POISSON, GBM_GAMMA, GBM_TWEEDIE)
         pred = exp(pred)
      case (GBM_PAIRWISE)
         do i = 1, size(pred)
            pred(i) = sigmoid(pred(i))
         end do
      case default
         continue
      end select
   end function gbm_predict_response

   function gbm_predict_trees(model, x) result(contribution)
      type(gbm_model), intent(in) :: model
      real(dp), intent(in) :: x(:, :)
      real(dp), allocatable :: contribution(:, :)
      integer :: i, t

      if (size(x, 2) /= model%n_features) error stop "gbm_predict_trees: feature count mismatch"
      allocate(contribution(size(x, 1), model%n_trees))
      do t = 1, model%n_trees
         do i = 1, size(x, 1)
            contribution(i, t) = predict_tree(model%trees(t), x(i, :))
         end do
      end do
   end function gbm_predict_trees

   function gbm_predict_staged(model, x, tree_counts, offset) result(pred)
      type(gbm_model), intent(in) :: model
      real(dp), intent(in) :: x(:, :)
      integer, intent(in) :: tree_counts(:)
      real(dp), intent(in), optional :: offset(:)
      real(dp), allocatable :: pred(:, :)
      real(dp), allocatable :: current(:)
      integer :: n, i, j, t, max_tree

      n = size(x, 1)
      if (size(x, 2) /= model%n_features) error stop "gbm_predict_staged: feature count mismatch"
      if (size(tree_counts) < 1) error stop "gbm_predict_staged: tree_counts cannot be empty"
      if (any(tree_counts < 0) .or. any(tree_counts > model%n_trees)) &
         error stop "gbm_predict_staged: invalid tree count"
      allocate(pred(n, size(tree_counts)), current(n))
      current = model%init_f
      if (present(offset)) then
         if (size(offset) /= n) error stop "gbm_predict_staged: offset length mismatch"
         current = current + offset
      end if
      do j = 1, size(tree_counts)
         if (tree_counts(j) == 0) pred(:, j) = current
      end do
      max_tree = maxval(tree_counts)
      do t = 1, max_tree
         do i = 1, n
            current(i) = current(i) + predict_tree(model%trees(t), x(i, :))
         end do
         do j = 1, size(tree_counts)
            if (tree_counts(j) == t) pred(:, j) = current
         end do
      end do
   end function gbm_predict_staged

   function gbm_predict_response_staged(model, x, tree_counts, offset) result(pred)
      type(gbm_model), intent(in) :: model
      real(dp), intent(in) :: x(:, :)
      integer, intent(in) :: tree_counts(:)
      real(dp), intent(in), optional :: offset(:)
      real(dp), allocatable :: pred(:, :)
      integer :: i, j

      pred = gbm_predict_staged(model, x, tree_counts, offset)
      select case (model%options%distribution)
      case (GBM_BERNOULLI)
         do j = 1, size(pred, 2)
            do i = 1, size(pred, 1)
               pred(i, j) = sigmoid(pred(i, j))
            end do
         end do
      case (GBM_ADABOOST)
         do j = 1, size(pred, 2)
            do i = 1, size(pred, 1)
               pred(i, j) = sigmoid(2.0_dp * pred(i, j))
            end do
         end do
      case (GBM_POISSON, GBM_GAMMA, GBM_TWEEDIE)
         pred = exp(pred)
      case (GBM_PAIRWISE)
         do j = 1, size(pred, 2)
            do i = 1, size(pred, 1)
               pred(i, j) = sigmoid(pred(i, j))
            end do
         end do
      case default
         continue
      end select
   end function gbm_predict_response_staged

   function gbm_partial_dependence(model, values, variables, n_trees) result(pred)
      type(gbm_model), intent(in) :: model
      real(dp), intent(in) :: values(:, :)
      integer, intent(in) :: variables(:)
      integer, intent(in), optional :: n_trees
      real(dp), allocatable :: pred(:)
      integer :: i, t, nt

      if (size(values, 2) /= size(variables)) error stop "gbm_partial_dependence: shape mismatch"
      if (size(variables) == 0) error stop "gbm_partial_dependence: variables cannot be empty"
      if (any(variables < 1) .or. any(variables > model%n_features)) &
         error stop "gbm_partial_dependence: variable index out of range"
      do i = 1, size(variables)
         if (count(variables == variables(i)) /= 1) error stop "gbm_partial_dependence: duplicate variable"
      end do
      nt = model%n_trees
      if (present(n_trees)) nt = n_trees
      if (nt < 0 .or. nt > model%n_trees) error stop "gbm_partial_dependence: invalid n_trees"

      allocate(pred(size(values, 1)))
      pred = model%init_f
      do t = 1, nt
         do i = 1, size(values, 1)
            pred(i) = pred(i) + partial_tree_prediction(model%trees(t), 1, values(i, :), variables)
         end do
      end do
   end function gbm_partial_dependence

   recursive real(dp) function partial_tree_prediction(tree, node_id, row, variables) result(pred)
      type(gbm_tree), intent(in) :: tree
      integer, intent(in) :: node_id
      real(dp), intent(in) :: row(:)
      integer, intent(in) :: variables(:)
      integer :: k, child, cat, l, r
      real(dp) :: xv, den

      if (tree%nodes(node_id)%is_terminal) then
         pred = tree%shrinkage * tree%nodes(node_id)%prediction
         return
      end if
      k = find_variable(tree%nodes(node_id)%split_var, variables)
      if (k > 0) then
         xv = row(k)
         if (ieee_is_nan(xv)) then
            child = tree%nodes(node_id)%missing
         else if (allocated(tree%nodes(node_id)%left_categories)) then
            cat = int(xv)
            if (any(tree%nodes(node_id)%left_categories == cat)) then
               child = tree%nodes(node_id)%left
            else
               child = tree%nodes(node_id)%right
            end if
         else if (xv < tree%nodes(node_id)%split_value) then
            child = tree%nodes(node_id)%left
         else
            child = tree%nodes(node_id)%right
         end if
         pred = partial_tree_prediction(tree, child, row, variables)
      else
         ! Match gbm_plot: marginalize an unselected split over left/right only,
         ! weighted by child training weights. The missing branch is not used.
         l = tree%nodes(node_id)%left
         r = tree%nodes(node_id)%right
         den = tree%nodes(l)%total_weight + tree%nodes(r)%total_weight
         if (den <= 0.0_dp) then
            pred = tree%shrinkage * tree%nodes(node_id)%prediction
         else
            pred = tree%nodes(l)%total_weight / den * partial_tree_prediction(tree, l, row, variables) + &
                   tree%nodes(r)%total_weight / den * partial_tree_prediction(tree, r, row, variables)
         end if
      end if
   end function partial_tree_prediction

   integer function find_variable(variable, variables) result(pos)
      integer, intent(in) :: variable, variables(:)
      integer :: i
      pos = 0
      do i = 1, size(variables)
         if (variables(i) == variable) then
            pos = i
            return
         end if
      end do
   end function find_variable

   subroutine gbm_relative_influence(model, influence, normalize, n_trees)
      type(gbm_model), intent(in) :: model
      real(dp), allocatable, intent(out) :: influence(:)
      logical, intent(in), optional :: normalize
      integer, intent(in), optional :: n_trees
      logical :: do_normalize
      integer :: t, nt
      real(dp) :: s

      nt = model%n_trees
      if (present(n_trees)) nt = n_trees
      if (nt < 0 .or. nt > model%n_trees) error stop "gbm_relative_influence: invalid n_trees"
      allocate(influence(model%n_features))
      influence = 0.0_dp
      do t = 1, nt
         call tree_relative_influence(model%trees(t), influence)
      end do
      do_normalize = .false.
      if (present(normalize)) do_normalize = normalize
      if (do_normalize) then
         s = sum(influence)
         if (s > 0.0_dp .and. ieee_is_finite(s)) influence = 100.0_dp * influence / s
      end if
   end subroutine gbm_relative_influence

   subroutine gbm_set_seed(seed)
      integer, intent(in) :: seed
      integer :: n, i
      integer, allocatable :: put(:)
      call random_seed(size=n)
      allocate(put(n))
      do i = 1, n
         put(i) = modulo(seed + 104729 * i, huge(1) - 1)
         if (put(i) == 0) put(i) = i
      end do
      call random_seed(put=put)
   end subroutine gbm_set_seed

end module gbm3_core
