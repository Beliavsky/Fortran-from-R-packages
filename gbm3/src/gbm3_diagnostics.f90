! SPDX-License-Identifier: GPL-2.0-or-later
! Modern Fortran translation of gbm3 3.0.3; see NOTICE.md and upstream/.
module gbm3_diagnostics
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
   use gbm3_kinds, only : dp
   use gbm3_constants
   use gbm3_types, only : gbm_model, gbm_tree, gbm_node
   use gbm3_math, only : shuffle_int, quiet_nan
   use gbm3_core, only : gbm_predict, gbm_partial_dependence
   use gbm3_distributions, only : dist_deviance
   use gbm3_pairwise, only : pairwise_deviance
   use gbm3_cox, only : cox_deviance
   implicit none
   private

   interface gbm_permutation_importance
      module procedure gbm_permutation_importance_vector
      module procedure gbm_permutation_importance_cox
   end interface gbm_permutation_importance

   public :: gbm_best_iteration
   public :: gbm_permutation_importance, gbm_permutation_importance_vector, gbm_permutation_importance_cox
   public :: gbm_get_tree, gbm_get_node, gbm_num_nodes
   public :: gbm_interaction_strength

   type :: pd_subset
      integer, allocatable :: variables(:)
      integer, allocatable :: positions(:)
      real(dp), allocatable :: values(:, :)
      real(dp), allocatable :: prediction(:)
      real(dp), allocatable :: count(:)
      integer :: sign = 1
   end type pd_subset

contains

   integer function gbm_best_iteration(model, method) result(iteration)
      type(gbm_model), intent(in) :: model
      character(*), intent(in) :: method
      real(dp) :: best, cumulative
      integer :: i
      character(len=len(method)) :: key

      iteration = 0
      best = huge(1.0_dp)
      key = lower_ascii(trim(adjustl(method)))
      select case (trim(key))
      case ('train')
         do i = 1, model%n_trees
            if (.not. ieee_is_finite(model%train_error(i))) cycle
            if (iteration == 0 .or. model%train_error(i) < best) then
               iteration = i
               best = model%train_error(i)
            end if
         end do
      case ('validation', 'valid', 'test')
         do i = 1, model%n_trees
            if (.not. ieee_is_finite(model%validation_error(i))) cycle
            if (iteration == 0 .or. model%validation_error(i) < best) then
               iteration = i
               best = model%validation_error(i)
            end if
         end do
      case ('oob_raw')
         ! Upstream gbm3 smooths OOB improvement with R's loess before taking
         ! the maximum cumulative improvement.  This explicitly named raw
         ! selector omits that R-only smoothing step.
         cumulative = 0.0_dp
         best = -huge(1.0_dp)
         do i = 1, model%n_trees
            if (.not. ieee_is_finite(model%oob_improvement(i))) cycle
            cumulative = cumulative + model%oob_improvement(i)
            if (iteration == 0 .or. cumulative > best) then
               iteration = i
               best = cumulative
            end if
         end do
      case default
         error stop "gbm_best_iteration: method must be train, validation/test, or oob_raw"
      end select
      if (iteration == 0) error stop "gbm_best_iteration: no finite performance values for requested method"
   end function gbm_best_iteration

   subroutine gbm_permutation_importance_vector(model, x, y, importance, n_trees, weight, offset, group, rescale)
      type(gbm_model), intent(in) :: model
      real(dp), intent(in) :: x(:, :), y(:)
      real(dp), allocatable, intent(out) :: importance(:)
      integer, intent(in), optional :: n_trees
      real(dp), intent(in), optional :: weight(:), offset(:)
      integer, intent(in), optional :: group(:)
      logical, intent(in), optional :: rescale

      real(dp), allocatable :: w(:), off(:), xperm(:, :), pred(:)
      integer, allocatable :: rows(:), grp(:)
      logical, allocatable :: used(:)
      integer :: n, p, nt, i, j
      real(dp) :: baseline, shuffled_loss, scale
      logical :: do_rescale

      n = size(x, 1)
      p = size(x, 2)
      if (size(y) /= n) error stop "gbm_permutation_importance: response length mismatch"
      if (p /= model%n_features) error stop "gbm_permutation_importance: feature count mismatch"
      if (model%options%distribution == GBM_COXPH) &
         error stop "gbm_permutation_importance: use the survival-response overload for Cox PH"
      nt = model%n_trees
      if (present(n_trees)) nt = n_trees
      if (nt < 1 .or. nt > model%n_trees) error stop "gbm_permutation_importance: invalid n_trees"

      allocate(w(n), off(n), grp(n), rows(n), xperm(n, p), used(p), importance(p))
      w = 1.0_dp
      if (present(weight)) then
         if (size(weight) /= n) error stop "gbm_permutation_importance: weight length mismatch"
         w = weight
      end if
      off = 0.0_dp
      if (present(offset)) then
         if (size(offset) /= n) error stop "gbm_permutation_importance: offset length mismatch"
         off = offset
      end if
      grp = [(i, i=1,n)]
      if (model%options%distribution == GBM_PAIRWISE) then
         if (.not. present(group)) error stop "gbm_permutation_importance: pairwise distribution requires group"
         if (size(group) /= n) error stop "gbm_permutation_importance: group length mismatch"
         grp = group
      end if

      pred = gbm_predict(model, x, n_trees=nt)
      if (model%options%distribution == GBM_PAIRWISE) then
         baseline = pairwise_deviance(y, grp, off, pred, w, model%options)
      else
         baseline = dist_deviance(y, off, w, pred, model%options)
      end if

      call mark_used_variables(model, nt, used)
      rows = [(i, i=1,n)]
      call shuffle_int(rows)
      xperm = x
      importance = 0.0_dp
      do j = 1, p
         if (.not. used(j)) cycle
         xperm(:, j) = x(rows, j)
         pred = gbm_predict(model, xperm, n_trees=nt)
         if (model%options%distribution == GBM_PAIRWISE) then
            shuffled_loss = pairwise_deviance(y, grp, off, pred, w, model%options)
         else
            shuffled_loss = dist_deviance(y, off, w, pred, model%options)
         end if
         importance(j) = shuffled_loss - baseline
         xperm(:, j) = x(:, j)
      end do

      do_rescale = .false.
      if (present(rescale)) do_rescale = rescale
      if (do_rescale) then
         scale = maxval(importance)
         if (abs(scale) > tiny(1.0_dp)) importance = importance / scale
      end if
   end subroutine gbm_permutation_importance_vector

   subroutine gbm_permutation_importance_cox(model, x, surv, importance, n_trees, weight, offset, strata, rescale)
      type(gbm_model), intent(in) :: model
      real(dp), intent(in) :: x(:, :), surv(:, :)
      real(dp), allocatable, intent(out) :: importance(:)
      integer, intent(in), optional :: n_trees
      real(dp), intent(in), optional :: weight(:), offset(:)
      integer, intent(in), optional :: strata(:)
      logical, intent(in), optional :: rescale

      real(dp), allocatable :: w(:), off(:), xperm(:, :), pred(:)
      integer, allocatable :: rows(:), str(:)
      logical, allocatable :: used(:)
      integer :: n, p, nt, i, j
      real(dp) :: baseline, shuffled_loss, scale
      logical :: do_rescale

      n = size(x, 1)
      p = size(x, 2)
      if (size(surv, 1) /= n) error stop "gbm_permutation_importance: survival response row mismatch"
      if (size(surv, 2) /= 2 .and. size(surv, 2) /= 3) &
         error stop "gbm_permutation_importance: survival response must have 2 or 3 columns"
      if (p /= model%n_features) error stop "gbm_permutation_importance: feature count mismatch"
      if (model%options%distribution /= GBM_COXPH) error stop "gbm_permutation_importance: model is not Cox PH"
      nt = model%n_trees
      if (present(n_trees)) nt = n_trees
      if (nt < 1 .or. nt > model%n_trees) error stop "gbm_permutation_importance: invalid n_trees"

      allocate(w(n), off(n), str(n), rows(n), xperm(n, p), used(p), importance(p))
      w = 1.0_dp
      if (present(weight)) then
         if (size(weight) /= n) error stop "gbm_permutation_importance: weight length mismatch"
         w = weight
      end if
      off = 0.0_dp
      if (present(offset)) then
         if (size(offset) /= n) error stop "gbm_permutation_importance: offset length mismatch"
         off = offset
      end if
      str = 1
      if (present(strata)) then
         if (size(strata) /= n) error stop "gbm_permutation_importance: strata length mismatch"
         str = strata
      end if

      pred = gbm_predict(model, x, n_trees=nt)
      baseline = cox_deviance(surv, str, off, pred, w, model%options)
      call mark_used_variables(model, nt, used)
      rows = [(i, i=1,n)]
      call shuffle_int(rows)
      xperm = x
      importance = 0.0_dp
      do j = 1, p
         if (.not. used(j)) cycle
         xperm(:, j) = x(rows, j)
         pred = gbm_predict(model, xperm, n_trees=nt)
         shuffled_loss = cox_deviance(surv, str, off, pred, w, model%options)
         importance(j) = shuffled_loss - baseline
         xperm(:, j) = x(:, j)
      end do

      do_rescale = .false.
      if (present(rescale)) do_rescale = rescale
      if (do_rescale) then
         scale = maxval(importance)
         if (abs(scale) > tiny(1.0_dp)) importance = importance / scale
      end if
   end subroutine gbm_permutation_importance_cox

   subroutine mark_used_variables(model, n_trees, used)
      type(gbm_model), intent(in) :: model
      integer, intent(in) :: n_trees
      logical, intent(out) :: used(model%n_features)
      integer :: t, node, variable
      used = .false.
      do t = 1, n_trees
         do node = 1, model%trees(t)%n_nodes
            if (model%trees(t)%nodes(node)%is_terminal) cycle
            variable = model%trees(t)%nodes(node)%split_var
            if (variable >= 1 .and. variable <= model%n_features) used(variable) = .true.
         end do
      end do
   end subroutine mark_used_variables

   real(dp) function gbm_interaction_strength(model, x, variables, n_trees) result(h)
      type(gbm_model), intent(in) :: model
      real(dp), intent(in) :: x(:, :)
      integer, intent(in) :: variables(:)
      integer, intent(in), optional :: n_trees

      type(pd_subset), allocatable :: subsets(:)
      logical, allocatable :: selected(:)
      real(dp), allocatable :: interaction(:), full_sq(:)
      integer, allocatable :: pos(:), vars(:)
      integer :: k, nmask, mask, full_mask, nt, i, j, row, match, nfull
      real(dp) :: numerator, denominator, mean_pred

      if (size(x, 2) /= model%n_features) error stop "gbm_interaction_strength: feature count mismatch"
      k = size(variables)
      if (k < 2) error stop "gbm_interaction_strength: at least two variables are required"
      if (k > model%options%interaction_depth) &
         error stop "gbm_interaction_strength: number of variables exceeds fitted interaction depth"
      if (k >= bit_size(1) - 1) error stop "gbm_interaction_strength: too many variables"
      if (any(variables < 1) .or. any(variables > model%n_features)) &
         error stop "gbm_interaction_strength: variable index out of range"
      do i = 1, k
         if (count(variables == variables(i)) /= 1) error stop "gbm_interaction_strength: duplicate variable"
      end do
      nt = model%n_trees
      if (present(n_trees)) nt = n_trees
      if (nt < 1 .or. nt > model%n_trees) error stop "gbm_interaction_strength: invalid n_trees"

      full_mask = 2 ** k - 1
      nmask = full_mask
      allocate(subsets(nmask), selected(k))
      do mask = 1, nmask
         do j = 1, k
            selected(j) = btest(mask, j - 1)
         end do
         pos = pack([(j, j=1,k)], selected)
         vars = pack(variables, selected)
         subsets(mask)%positions = pos
         subsets(mask)%variables = vars
         call unique_rows_with_counts(x(:, vars), subsets(mask)%values, subsets(mask)%count)
         subsets(mask)%prediction = gbm_partial_dependence(model, subsets(mask)%values, vars, nt)
         mean_pred = sum(subsets(mask)%count * subsets(mask)%prediction) / sum(subsets(mask)%count)
         subsets(mask)%prediction = subsets(mask)%prediction - mean_pred
         if (mod(size(vars), 2) == mod(k, 2)) then
            subsets(mask)%sign = 1
         else
            subsets(mask)%sign = -1
         end if
      end do

      nfull = size(subsets(full_mask)%prediction)
      allocate(interaction(nfull), full_sq(nfull))
      interaction = subsets(full_mask)%prediction
      full_sq = subsets(full_mask)%prediction ** 2
      do mask = 1, full_mask - 1
         do row = 1, nfull
            match = find_matching_projection(subsets(mask)%values, subsets(full_mask)%values(row, :), &
                                             subsets(mask)%positions)
            if (match == 0) error stop "gbm_interaction_strength: internal projection match failure"
            interaction(row) = interaction(row) + real(subsets(mask)%sign, dp) * subsets(mask)%prediction(match)
         end do
      end do

      numerator = sum(subsets(full_mask)%count * interaction ** 2) / sum(subsets(full_mask)%count)
      denominator = sum(subsets(full_mask)%count * full_sq) / sum(subsets(full_mask)%count)
      if (denominator <= tiny(1.0_dp) .or. numerator < 0.0_dp) then
         h = quiet_nan()
      else if (numerator / denominator > 1.0_dp) then
         h = quiet_nan()
      else
         h = sqrt(max(0.0_dp, numerator / denominator))
      end if
   end function gbm_interaction_strength

   subroutine unique_rows_with_counts(x, unique_x, counts)
      real(dp), intent(in) :: x(:, :)
      real(dp), allocatable, intent(out) :: unique_x(:, :), counts(:)
      real(dp), allocatable :: work(:, :), work_count(:)
      integer :: i, j, n_unique
      logical :: found

      allocate(work(size(x, 1), size(x, 2)), work_count(size(x, 1)))
      n_unique = 0
      work_count = 0.0_dp
      do i = 1, size(x, 1)
         found = .false.
         do j = 1, n_unique
            if (rows_identical(x(i, :), work(j, :))) then
               work_count(j) = work_count(j) + 1.0_dp
               found = .true.
               exit
            end if
         end do
         if (.not. found) then
            n_unique = n_unique + 1
            work(n_unique, :) = x(i, :)
            work_count(n_unique) = 1.0_dp
         end if
      end do
      allocate(unique_x(n_unique, size(x, 2)), counts(n_unique))
      unique_x = work(1:n_unique, :)
      counts = work_count(1:n_unique)
   end subroutine unique_rows_with_counts

   integer function find_matching_projection(values, full_row, positions) result(idx)
      real(dp), intent(in) :: values(:, :), full_row(:)
      integer, intent(in) :: positions(:)
      integer :: i
      idx = 0
      do i = 1, size(values, 1)
         if (rows_identical(values(i, :), full_row(positions))) then
            idx = i
            return
         end if
      end do
   end function find_matching_projection

   logical function rows_identical(a, b) result(same)
      real(dp), intent(in) :: a(:), b(:)
      integer :: i
      if (size(a) /= size(b)) then
         same = .false.
         return
      end if
      same = .true.
      do i = 1, size(a)
         if (ieee_is_nan(a(i)) .and. ieee_is_nan(b(i))) cycle
         if (ieee_is_nan(a(i)) .or. ieee_is_nan(b(i))) then
            same = .false.
            return
         end if
         if (.not. (a(i) <= b(i) .and. a(i) >= b(i))) then
            same = .false.
            return
         end if
      end do
   end function rows_identical

   pure function lower_ascii(text) result(lower)
      character(*), intent(in) :: text
      character(len=len(text)) :: lower
      integer :: i, code
      lower = text
      do i = 1, len(text)
         code = iachar(text(i:i))
         if (code >= iachar('A') .and. code <= iachar('Z')) lower(i:i) = achar(code + iachar('a') - iachar('A'))
      end do
   end function lower_ascii

   integer function gbm_num_nodes(model, tree_index) result(n_nodes)
      type(gbm_model), intent(in) :: model
      integer, intent(in) :: tree_index
      call validate_tree_index(model, tree_index)
      n_nodes = model%trees(tree_index)%n_nodes
   end function gbm_num_nodes

   function gbm_get_tree(model, tree_index) result(tree)
      type(gbm_model), intent(in) :: model
      integer, intent(in) :: tree_index
      type(gbm_tree) :: tree
      call validate_tree_index(model, tree_index)
      tree = model%trees(tree_index)
   end function gbm_get_tree

   function gbm_get_node(model, tree_index, node_index) result(node)
      type(gbm_model), intent(in) :: model
      integer, intent(in) :: tree_index, node_index
      type(gbm_node) :: node
      call validate_tree_index(model, tree_index)
      if (node_index < 1 .or. node_index > model%trees(tree_index)%n_nodes) &
         error stop "gbm_get_node: node index out of range"
      node = model%trees(tree_index)%nodes(node_index)
   end function gbm_get_node

   subroutine validate_tree_index(model, tree_index)
      type(gbm_model), intent(in) :: model
      integer, intent(in) :: tree_index
      if (tree_index < 1 .or. tree_index > model%n_trees) error stop "gbm tree index out of range"
   end subroutine validate_tree_index

end module gbm3_diagnostics
