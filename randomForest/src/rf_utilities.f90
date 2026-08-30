! SPDX-License-Identifier: GPL-2.0-or-later
module rf_utilities
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
   use r_kinds, only : dp
   use r_quantiles, only : r_median
   use r_robust, only : r_mad
   use rf_types, only : rf_classification_forest, rf_regression_forest, RF_INTERIOR, RF_TERMINAL
   use rf_classification, only : predict_classification
   use rf_regression, only : predict_regression
   implicit none
   private

   interface tree_sizes
      module procedure class_tree_sizes
      module procedure reg_tree_sizes
   end interface tree_sizes

   interface variable_usage
      module procedure class_variable_usage
      module procedure reg_variable_usage
   end interface variable_usage

   public :: classification_margin, outlier_scores, roughfix_numeric
   public :: class_centers, tree_sizes, variable_usage
   public :: partial_dependence_regression, partial_dependence_classification

contains

   subroutine classification_margin(votes, observed, margin)
      real(dp), intent(in) :: votes(:,:)
      integer, intent(in) :: observed(:)
      real(dp), intent(out) :: margin(:)
      integer :: i, c
      real(dp) :: total, correct, other

      if (size(votes, 1) /= size(observed) .or. size(margin) /= size(observed)) &
         error stop 'classification_margin: incompatible shapes'
      do i = 1, size(observed)
         if (observed(i) < 1 .or. observed(i) > size(votes, 2)) error stop 'classification_margin: invalid class label'
         total = sum(votes(i, :))
         if (total <= 0.0_dp) then
            margin(i) = 0.0_dp
            cycle
         end if
         correct = votes(i, observed(i)) / total
         other = -huge(1.0_dp)
         do c = 1, size(votes, 2)
            if (c /= observed(i)) other = max(other, votes(i, c) / total)
         end do
         margin(i) = correct - other
      end do
   end subroutine classification_margin

   subroutine outlier_scores(proximity, class_label, score)
      real(dp), intent(in) :: proximity(:,:)
      integer, intent(in) :: class_label(:)
      real(dp), intent(out) :: score(:)
      real(dp), allocatable :: raw(:), class_raw(:)
      integer :: n, c, i, j, m, pos
      real(dp) :: s, med, mad

      n = size(proximity, 1)
      if (size(proximity, 2) /= n .or. size(class_label) /= n .or. size(score) /= n) &
         error stop 'outlier_scores: incompatible shapes'
      allocate(raw(n))
      raw = 0.0_dp
      do c = 1, maxval(class_label)
         m = count(class_label == c)
         if (m <= 0) cycle
         do i = 1, n
            if (class_label(i) /= c) cycle
            s = 0.0_dp
            do j = 1, n
               if (class_label(j) == c) s = s + proximity(i, j) ** 2
            end do
            if (s <= 0.0_dp) s = 1.0_dp
            raw(i) = real(n, dp) / s
         end do
         allocate(class_raw(m))
         pos = 0
         do i = 1, n
            if (class_label(i) == c) then
               pos = pos + 1
               class_raw(pos) = raw(i)
            end if
         end do
         med = r_median(class_raw)
         mad = r_mad(class_raw)
         if (mad > 0.0_dp) then
            do i = 1, n
               if (class_label(i) == c) score(i) = (raw(i) - med) / mad
            end do
         else
            do i = 1, n
               if (class_label(i) == c) score(i) = 0.0_dp
            end do
         end if
         deallocate(class_raw)
      end do
   end subroutine outlier_scores

   subroutine roughfix_numeric(x)
      real(dp), intent(inout) :: x(:,:)
      integer :: j
      real(dp) :: med
      logical, allocatable :: missing(:)

      allocate(missing(size(x, 1)))
      do j = 1, size(x, 2)
         missing = ieee_is_nan(x(:, j))
         if (any(missing)) then
            med = r_median(x(:, j), na_rm=.true.)
            where (missing) x(:, j) = med
         end if
      end do
   end subroutine roughfix_numeric

   subroutine class_centers(x, label, proximity, centers, n_neighbors)
      real(dp), intent(in) :: x(:,:), proximity(:,:)
      integer, intent(in) :: label(:)
      real(dp), intent(out) :: centers(:,:)
      integer, intent(in), optional :: n_neighbors
      integer, allocatable :: order(:), neighbor_label_count(:,:), class_mode(:), selected(:)
      real(dp), allocatable :: values(:)
      integer :: n, p, nclass, nn, i, j, c, k, best_case, best_count, nsel

      n = size(x, 1)
      p = size(x, 2)
      nclass = maxval(label)
      if (size(label) /= n .or. size(proximity, 1) /= n .or. size(proximity, 2) /= n) &
         error stop 'class_centers: incompatible shapes'
      if (size(centers, 1) /= nclass .or. size(centers, 2) /= p) error stop 'class_centers: centers has wrong shape'
      if (present(n_neighbors)) then
         nn = n_neighbors
      else
         nn = minval([(count(label == c), c = 1, nclass)]) - 1
      end if
      nn = max(1, min(nn, n))
      allocate(order(n), neighbor_label_count(n, nclass), class_mode(nclass), selected(n), values(n))
      neighbor_label_count = 0

      do i = 1, n
         order = [(j, j = 1, n)]
         call sort_neighbors_desc(proximity(i, :), order)
         do k = 1, nn
            c = label(order(k))
            neighbor_label_count(i, c) = neighbor_label_count(i, c) + 1
         end do
      end do

      do c = 1, nclass
         best_case = 1
         best_count = -1
         do i = 1, n
            if (neighbor_label_count(i, c) > best_count) then
               best_count = neighbor_label_count(i, c)
               best_case = i
            end if
         end do
         class_mode(c) = best_case
         order = [(j, j = 1, n)]
         call sort_neighbors_desc(proximity(best_case, :), order)
         nsel = 0
         do k = 1, nn
            if (label(order(k)) == c) then
               nsel = nsel + 1
               selected(nsel) = order(k)
            end if
         end do
         if (nsel == 0) then
            centers(c, :) = x(best_case, :)
         else
            do j = 1, p
               values(1:nsel) = x(selected(1:nsel), j)
               centers(c, j) = r_median(values(1:nsel))
            end do
         end if
      end do
   end subroutine class_centers

   subroutine sort_neighbors_desc(values, order)
      real(dp), intent(in) :: values(:)
      integer, intent(inout) :: order(:)
      integer :: i, j, key
      real(dp) :: key_value

      do i = 2, size(order)
         key = order(i)
         key_value = values(key)
         j = i - 1
         do while (j >= 1)
            if (values(order(j)) >= key_value) exit
            order(j + 1) = order(j)
            j = j - 1
         end do
         order(j + 1) = key
      end do
   end subroutine sort_neighbors_desc

   subroutine class_tree_sizes(forest, sizes, terminal_only)
      type(rf_classification_forest), intent(in) :: forest
      integer, intent(out) :: sizes(:)
      logical, intent(in), optional :: terminal_only
      logical :: terminal
      integer :: t

      if (size(sizes) /= size(forest%trees)) error stop 'tree_sizes: sizes has wrong length'
      terminal = .true.
      if (present(terminal_only)) terminal = terminal_only
      do t = 1, size(forest%trees)
         if (terminal) then
            sizes(t) = count(forest%trees(t)%status(1:forest%trees(t)%n_nodes) == RF_TERMINAL)
         else
            sizes(t) = forest%trees(t)%n_nodes
         end if
      end do
   end subroutine class_tree_sizes

   subroutine reg_tree_sizes(forest, sizes, terminal_only)
      type(rf_regression_forest), intent(in) :: forest
      integer, intent(out) :: sizes(:)
      logical, intent(in), optional :: terminal_only
      logical :: terminal
      integer :: t

      if (size(sizes) /= size(forest%trees)) error stop 'tree_sizes: sizes has wrong length'
      terminal = .true.
      if (present(terminal_only)) terminal = terminal_only
      do t = 1, size(forest%trees)
         if (terminal) then
            sizes(t) = count(forest%trees(t)%status(1:forest%trees(t)%n_nodes) == RF_TERMINAL)
         else
            sizes(t) = forest%trees(t)%n_nodes
         end if
      end do
   end subroutine reg_tree_sizes

   subroutine class_variable_usage(forest, counts, by_tree)
      type(rf_classification_forest), intent(in) :: forest
      integer, intent(out) :: counts(:,:)
      logical, intent(in), optional :: by_tree
      logical :: separate
      integer :: t, node, v

      separate = .true.
      if (present(by_tree)) separate = by_tree
      counts = 0
      if (separate) then
         if (size(counts, 1) /= forest%nvar .or. size(counts, 2) /= size(forest%trees)) &
            error stop 'variable_usage: counts must be nvar by ntree'
         do t = 1, size(forest%trees)
            do node = 1, forest%trees(t)%n_nodes
               if (forest%trees(t)%status(node) == RF_INTERIOR) then
                  v = forest%trees(t)%split_var(node)
                  counts(v, t) = counts(v, t) + 1
               end if
            end do
         end do
      else
         if (size(counts, 1) /= forest%nvar .or. size(counts, 2) /= 1) &
            error stop 'variable_usage: aggregate counts must be nvar by 1'
         do t = 1, size(forest%trees)
            do node = 1, forest%trees(t)%n_nodes
               if (forest%trees(t)%status(node) == RF_INTERIOR) then
                  v = forest%trees(t)%split_var(node)
                  counts(v, 1) = counts(v, 1) + 1
               end if
            end do
         end do
      end if
   end subroutine class_variable_usage

   subroutine reg_variable_usage(forest, counts, by_tree)
      type(rf_regression_forest), intent(in) :: forest
      integer, intent(out) :: counts(:,:)
      logical, intent(in), optional :: by_tree
      logical :: separate
      integer :: t, node, v

      separate = .true.
      if (present(by_tree)) separate = by_tree
      counts = 0
      if (separate) then
         if (size(counts, 1) /= forest%nvar .or. size(counts, 2) /= size(forest%trees)) &
            error stop 'variable_usage: counts must be nvar by ntree'
         do t = 1, size(forest%trees)
            do node = 1, forest%trees(t)%n_nodes
               if (forest%trees(t)%status(node) == RF_INTERIOR) then
                  v = forest%trees(t)%split_var(node)
                  counts(v, t) = counts(v, t) + 1
               end if
            end do
         end do
      else
         if (size(counts, 1) /= forest%nvar .or. size(counts, 2) /= 1) &
            error stop 'variable_usage: aggregate counts must be nvar by 1'
         do t = 1, size(forest%trees)
            do node = 1, forest%trees(t)%n_nodes
               if (forest%trees(t)%status(node) == RF_INTERIOR) then
                  v = forest%trees(t)%split_var(node)
                  counts(v, 1) = counts(v, 1) + 1
               end if
            end do
         end do
      end if
   end subroutine reg_variable_usage

   subroutine partial_dependence_regression(forest, data, variable, grid, effect, weights)
      type(rf_regression_forest), intent(in) :: forest
      real(dp), intent(in) :: data(:,:), grid(:)
      integer, intent(in) :: variable
      real(dp), intent(out) :: effect(:)
      real(dp), intent(in), optional :: weights(:)
      real(dp), allocatable :: modified(:,:), pred(:), w(:)
      integer :: g

      if (size(effect) /= size(grid)) error stop 'partial_dependence_regression: effect has wrong size'
      if (variable < 1 .or. variable > size(data, 2)) error stop 'partial_dependence_regression: invalid variable'
      allocate(modified(size(data, 1), size(data, 2)), pred(size(data, 1)), w(size(data, 1)))
      if (present(weights)) then
         if (size(weights) /= size(data, 1) .or. any(weights < 0.0_dp) .or. sum(weights) <= 0.0_dp) &
            error stop 'partial_dependence_regression: invalid weights'
         w = weights / sum(weights)
      else
         w = 1.0_dp / real(size(data, 1), dp)
      end if
      do g = 1, size(grid)
         modified = data
         modified(:, variable) = grid(g)
         call predict_regression(forest, modified, pred)
         effect(g) = sum(w * pred)
      end do
   end subroutine partial_dependence_regression

   subroutine partial_dependence_classification(forest, data, variable, grid, focus_class, effect, weights)
      type(rf_classification_forest), intent(in) :: forest
      real(dp), intent(in) :: data(:,:), grid(:)
      integer, intent(in) :: variable, focus_class
      real(dp), intent(out) :: effect(:)
      real(dp), intent(in), optional :: weights(:)
      real(dp), allocatable :: modified(:,:), prob(:,:), w(:), log_prob(:,:)
      integer, allocatable :: pred(:)
      integer :: g, i
      real(dp) :: eps_prob

      if (size(effect) /= size(grid)) error stop 'partial_dependence_classification: effect has wrong size'
      if (variable < 1 .or. variable > size(data, 2)) error stop 'partial_dependence_classification: invalid variable'
      if (focus_class < 1 .or. focus_class > forest%nclass) error stop 'partial_dependence_classification: invalid class'
      allocate(modified(size(data, 1), size(data, 2)), prob(size(data, 1), forest%nclass))
      allocate(log_prob(size(data, 1), forest%nclass), pred(size(data, 1)), w(size(data, 1)))
      if (present(weights)) then
         if (size(weights) /= size(data, 1) .or. any(weights < 0.0_dp) .or. sum(weights) <= 0.0_dp) &
            error stop 'partial_dependence_classification: invalid weights'
         w = weights / sum(weights)
      else
         w = 1.0_dp / real(size(data, 1), dp)
      end if
      eps_prob = epsilon(1.0_dp)
      do g = 1, size(grid)
         modified = data
         modified(:, variable) = grid(g)
         call predict_classification(forest, modified, pred, probabilities=prob)
         log_prob = log(max(prob, eps_prob))
         effect(g) = 0.0_dp
         do i = 1, size(data, 1)
            effect(g) = effect(g) + w(i) * (log_prob(i, focus_class) - sum(log_prob(i, :)) / real(forest%nclass, dp))
         end do
      end do
   end subroutine partial_dependence_classification

end module rf_utilities
