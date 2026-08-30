! SPDX-License-Identifier: GPL-2.0-or-later
! Modern Fortran translation of gbm3 3.0.3; see NOTICE.md and upstream/.
module gbm3_tree
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
   use gbm3_kinds, only : dp
   use gbm3_math, only : argsort_real, shuffle_int
   use gbm3_types, only : gbm_tree, gbm_node, split_candidate, node_stat
   implicit none
   private
   public :: grow_tree, adjust_tree, predict_tree, tree_relative_influence

contains

   subroutine grow_tree(x, residual, weight, in_bag, var_classes, monotone, depth, min_obs, &
                        num_features, shrinkage, tree, assignment)
      real(dp), intent(in) :: x(:, :), residual(:), weight(:)
      logical, intent(in) :: in_bag(:)
      integer, intent(in) :: var_classes(:), monotone(:)
      integer, intent(in) :: depth, min_obs, num_features
      real(dp), intent(in) :: shrinkage
      type(gbm_tree), intent(out) :: tree
      integer, allocatable, intent(out) :: assignment(:)

      integer :: n, p, max_nodes, max_leaves, d, j, k, nfeat
      integer, allocatable :: terminals(:), feature_order(:)
      integer :: nterm, node_id, best_pos
      type(split_candidate), allocatable :: best(:)
      type(split_candidate) :: cand, chosen
      type(node_stat) :: rootstat
      real(dp) :: best_imp

      n = size(x, 1)
      p = size(x, 2)
      if (size(residual) /= n .or. size(weight) /= n .or. size(in_bag) /= n) &
         error stop "grow_tree: row shape mismatch"
      if (size(var_classes) /= p .or. size(monotone) /= p) &
         error stop "grow_tree: feature shape mismatch"
      if (depth < 1) error stop "grow_tree: interaction depth must be positive"

      max_nodes = 1 + 3 * depth
      max_leaves = 1 + 2 * depth
      allocate(tree%nodes(max_nodes), assignment(n), terminals(max_leaves))
      tree%n_nodes = 1
      tree%shrinkage = shrinkage
      assignment = 1
      terminals = 0
      terminals(1) = 1
      nterm = 1

      rootstat = node_stat()
      do j = 1, n
         if (in_bag(j)) then
            rootstat%n = rootstat%n + 1
            rootstat%wr = rootstat%wr + weight(j) * residual(j)
            rootstat%w = rootstat%w + weight(j)
         end if
      end do
      tree%nodes(1)%prediction = rootstat%prediction()
      tree%nodes(1)%total_weight = rootstat%w
      tree%nodes(1)%num_obs = rootstat%n

      allocate(feature_order(p))
      do j = 1, p
         feature_order(j) = j
      end do
      nfeat = num_features
      if (nfeat <= 0 .or. nfeat > p) nfeat = p

      do d = 1, depth
         call shuffle_int(feature_order)
         allocate(best(nterm))
         do j = 1, nterm
            best(j) = split_candidate()
            best(j)%node_id = terminals(j)
         end do

         ! Match gbm3's feature-order tie behavior: the first feature in the
         ! shuffled subset wins an exact improvement tie.
         do k = 1, nfeat
            j = feature_order(k)
            do best_pos = 1, nterm
               node_id = terminals(best_pos)
               if (var_classes(j) <= 0) then
                  call continuous_split(x(:, j), residual, weight, in_bag, assignment, node_id, &
                                        min_obs, monotone(j), cand)
               else
                  call categorical_split(x(:, j), residual, weight, in_bag, assignment, node_id, &
                                         min_obs, var_classes(j), cand)
               end if
               cand%split_var = j
               cand%split_class = max(0, var_classes(j))
               cand%node_id = node_id
               cand%bias = k
               if (cand%valid) then
                  if ((.not. best(best_pos)%valid) .or. cand%improvement > best(best_pos)%improvement) then
                     best(best_pos) = cand
                  end if
               end if
            end do
         end do

         best_imp = -huge(1.0_dp)
         best_pos = 0
         do j = 1, nterm
            if (best(j)%valid .and. best(j)%improvement > best_imp) then
               best_imp = best(j)%improvement
               best_pos = j
            end if
         end do
         if (best_pos == 0 .or. best_imp <= 0.0_dp) then
            deallocate(best)
            exit
         end if

         chosen = best(best_pos)
         node_id = terminals(best_pos)
         call split_node(tree, node_id, chosen)
         call reassign_rows(x, assignment, node_id, tree)

         ! Left child replaces the split node's terminal slot. Right and missing
         ! are appended, exactly as in the C++ terminal-node vector.
         terminals(best_pos) = tree%nodes(node_id)%left
         nterm = nterm + 1
         terminals(nterm) = tree%nodes(node_id)%right
         nterm = nterm + 1
         terminals(nterm) = tree%nodes(node_id)%missing

         deallocate(best)
      end do
   end subroutine grow_tree

   subroutine continuous_split(xcol, residual, weight, in_bag, assignment, node_id, &
                               min_obs, monotone, cand)
      real(dp), intent(in) :: xcol(:), residual(:), weight(:)
      logical, intent(in) :: in_bag(:)
      integer, intent(in) :: assignment(:), node_id, min_obs, monotone
      type(split_candidate), intent(out) :: cand

      integer, allocatable :: rows(:), ord(:)
      real(dp), allocatable :: vals(:)
      integer :: i, m, nr, row, nextrow
      type(node_stat) :: total, left, right, missing
      real(dp) :: last_x, xval, imp, parent_pred

      cand = split_candidate()
      nr = count(in_bag .and. assignment == node_id)
      if (nr <= 0) return
      allocate(rows(nr))
      m = 0
      total = node_stat()
      missing = node_stat()
      do i = 1, size(xcol)
         if (in_bag(i) .and. assignment(i) == node_id) then
            total%n = total%n + 1
            total%wr = total%wr + weight(i) * residual(i)
            total%w = total%w + weight(i)
            if (ieee_is_nan(xcol(i))) then
               missing%n = missing%n + 1
               missing%wr = missing%wr + weight(i) * residual(i)
               missing%w = missing%w + weight(i)
            else
               m = m + 1
               rows(m) = i
            end if
         end if
      end do
      if (m < 2 * min_obs) return

      allocate(vals(m), ord(m))
      do i = 1, m
         vals(i) = xcol(rows(i))
      end do
      call argsort_real(vals, ord)

      left = node_stat()
      right%n = total%n - missing%n
      right%wr = total%wr - missing%wr
      right%w = total%w - missing%w
      parent_pred = total%prediction()

      last_x = vals(ord(1))
      do i = 1, m - 1
         row = rows(ord(i))
         left%n = left%n + 1
         left%wr = left%wr + weight(row) * residual(row)
         left%w = left%w + weight(row)
         right%n = right%n - 1
         right%wr = right%wr - weight(row) * residual(row)
         right%w = right%w - weight(row)

         nextrow = rows(ord(i + 1))
         xval = xcol(nextrow)
         if ((last_x < xval .or. last_x > xval) .and. left%n >= min_obs .and. right%n >= min_obs) then
            if (monotone == 0 .or. real(monotone, dp) * unweighted_gradient(right, left) > 0.0_dp) then
               imp = node_improvement(left, right, missing)
               if ((.not. cand%valid) .or. imp > cand%improvement) then
                  cand%valid = .true.
                  cand%split_value = 0.5_dp * (last_x + xval)
                  cand%improvement = imp
                  cand%left = left
                  cand%right = right
                  cand%missing = missing
               end if
            end if
         end if
         last_x = xval
      end do

      if (cand%valid .and. cand%missing%n == 0) then
         cand%missing%n = 0
         cand%missing%w = total%w
         cand%missing%wr = parent_pred * total%w
      end if
   end subroutine continuous_split

   subroutine categorical_split(xcol, residual, weight, in_bag, assignment, node_id, &
                                min_obs, nclass, cand)
      real(dp), intent(in) :: xcol(:), residual(:), weight(:)
      logical, intent(in) :: in_bag(:)
      integer, intent(in) :: assignment(:), node_id, min_obs, nclass
      type(split_candidate), intent(out) :: cand

      type(node_stat), allocatable :: groups(:)
      type(node_stat) :: total, left, right, missing
      real(dp), allocatable :: means(:)
      integer, allocatable :: ord(:)
      integer :: i, cat, finite_count, row
      real(dp) :: imp, parent_pred

      cand = split_candidate()
      if (nclass < 2) return
      allocate(groups(nclass), means(nclass), ord(nclass))
      total = node_stat()
      missing = node_stat()
      groups = node_stat()

      do row = 1, size(xcol)
         if (.not. in_bag(row) .or. assignment(row) /= node_id) cycle
         total%n = total%n + 1
         total%wr = total%wr + weight(row) * residual(row)
         total%w = total%w + weight(row)
         if (ieee_is_nan(xcol(row))) then
            missing%n = missing%n + 1
            missing%wr = missing%wr + weight(row) * residual(row)
            missing%w = missing%w + weight(row)
         else
            cat = int(xcol(row)) + 1
            if (cat < 1 .or. cat > nclass) cycle
            groups(cat)%n = groups(cat)%n + 1
            groups(cat)%wr = groups(cat)%wr + weight(row) * residual(row)
            groups(cat)%w = groups(cat)%w + weight(row)
         end if
      end do

      finite_count = 0
      do i = 1, nclass
         if (groups(i)%w > 0.0_dp) then
            means(i) = groups(i)%prediction()
            finite_count = finite_count + 1
         else
            means(i) = huge(1.0_dp)
         end if
      end do
      if (finite_count < 2) return
      call argsort_real(means, ord)

      left = node_stat()
      right%n = total%n - missing%n
      right%wr = total%wr - missing%wr
      right%w = total%w - missing%w
      do i = 1, finite_count - 1
         cat = ord(i)
         left%n = left%n + groups(cat)%n
         left%wr = left%wr + groups(cat)%wr
         left%w = left%w + groups(cat)%w
         right%n = right%n - groups(cat)%n
         right%wr = right%wr - groups(cat)%wr
         right%w = right%w - groups(cat)%w
         if (left%n >= min_obs .and. right%n >= min_obs) then
            imp = node_improvement(left, right, missing)
            if ((.not. cand%valid) .or. imp > cand%improvement) then
               cand%valid = .true.
               cand%split_value = real(i - 1, dp)
               cand%improvement = imp
               cand%left = left
               cand%right = right
               cand%missing = missing
               if (allocated(cand%ordering)) deallocate(cand%ordering)
               allocate(cand%ordering(nclass))
               cand%ordering = ord - 1
            end if
         end if
      end do

      if (cand%valid .and. cand%missing%n == 0) then
         parent_pred = total%prediction()
         cand%missing%n = 0
         cand%missing%w = total%w
         cand%missing%wr = parent_pred * total%w
      end if
   end subroutine categorical_split

   pure real(dp) function unweighted_gradient(a, b) result(v)
      type(node_stat), intent(in) :: a, b
      v = a%wr * b%w - b%wr * a%w
   end function unweighted_gradient

   pure real(dp) function variance_reduction(a, b) result(v)
      type(node_stat), intent(in) :: a, b
      real(dp) :: d
      if (a%w <= 0.0_dp .or. b%w <= 0.0_dp) then
         v = 0.0_dp
      else
         d = a%prediction() - b%prediction()
         v = a%w * b%w * d * d
      end if
   end function variance_reduction

   pure real(dp) function node_improvement(left, right, missing) result(v)
      type(node_stat), intent(in) :: left, right, missing
      real(dp) :: den
      if (missing%n == 0) then
         den = left%w + right%w
         if (den <= 0.0_dp) then
            v = 0.0_dp
         else
            v = variance_reduction(left, right) / den
         end if
      else
         den = left%w + right%w + missing%w
         if (den <= 0.0_dp) then
            v = 0.0_dp
         else
            v = (variance_reduction(left, right) + variance_reduction(left, missing) + &
                 variance_reduction(right, missing)) / den
         end if
      end if
   end function node_improvement

   subroutine split_node(tree, node_id, cand)
      type(gbm_tree), intent(inout) :: tree
      integer, intent(in) :: node_id
      type(split_candidate), intent(in) :: cand
      integer :: l, r, m, nleft

      l = tree%n_nodes + 1
      r = tree%n_nodes + 2
      m = tree%n_nodes + 3
      tree%n_nodes = tree%n_nodes + 3

      tree%nodes(node_id)%is_terminal = .false.
      tree%nodes(node_id)%split_var = cand%split_var
      tree%nodes(node_id)%split_value = cand%split_value
      tree%nodes(node_id)%improvement = cand%improvement
      tree%nodes(node_id)%left = l
      tree%nodes(node_id)%right = r
      tree%nodes(node_id)%missing = m
      if (cand%split_class > 0) then
         nleft = 1 + int(cand%split_value)
         allocate(tree%nodes(node_id)%left_categories(nleft))
         tree%nodes(node_id)%left_categories = cand%ordering(1:nleft)
      end if

      call set_node_from_stat(tree%nodes(l), cand%left)
      call set_node_from_stat(tree%nodes(r), cand%right)
      call set_node_from_stat(tree%nodes(m), cand%missing)
   end subroutine split_node

   subroutine set_node_from_stat(node, stat)
      type(gbm_node), intent(inout) :: node
      type(node_stat), intent(in) :: stat
      node%is_terminal = .true.
      node%prediction = stat%prediction()
      node%total_weight = stat%w
      node%num_obs = stat%n
   end subroutine set_node_from_stat

   subroutine reassign_rows(x, assignment, split_node_id, tree)
      real(dp), intent(in) :: x(:, :)
      integer, intent(inout) :: assignment(:)
      integer, intent(in) :: split_node_id
      type(gbm_tree), intent(in) :: tree
      integer :: i, child
      do i = 1, size(assignment)
         if (assignment(i) == split_node_id) then
            child = choose_child(tree%nodes(split_node_id), x(i, :))
            assignment(i) = child
         end if
      end do
   end subroutine reassign_rows

   integer function choose_child(node, row) result(child)
      type(gbm_node), intent(in) :: node
      real(dp), intent(in) :: row(:)
      real(dp) :: xv
      integer :: cat
      xv = row(node%split_var)
      if (ieee_is_nan(xv)) then
         child = node%missing
      else if (allocated(node%left_categories)) then
         cat = int(xv)
         if (any(node%left_categories == cat)) then
            child = node%left
         else
            child = node%right
         end if
      else if (xv < node%split_value) then
         child = node%left
      else
         child = node%right
      end if
   end function choose_child

   recursive subroutine adjust_node(tree, idx, min_obs)
      type(gbm_tree), intent(inout) :: tree
      integer, intent(in) :: idx, min_obs
      integer :: l, r, m
      real(dp) :: den
      if (tree%nodes(idx)%is_terminal) return
      l = tree%nodes(idx)%left
      r = tree%nodes(idx)%right
      m = tree%nodes(idx)%missing
      call adjust_node(tree, l, min_obs)
      call adjust_node(tree, r, min_obs)

      if (tree%nodes(m)%is_terminal .and. tree%nodes(m)%num_obs < min_obs) then
         den = tree%nodes(l)%total_weight + tree%nodes(r)%total_weight
         if (den > 0.0_dp) then
            tree%nodes(idx)%prediction = (tree%nodes(l)%total_weight * tree%nodes(l)%prediction + &
                                          tree%nodes(r)%total_weight * tree%nodes(r)%prediction) / den
         end if
         tree%nodes(m)%prediction = tree%nodes(idx)%prediction
      else
         call adjust_node(tree, m, min_obs)
         den = tree%nodes(l)%total_weight + tree%nodes(r)%total_weight + tree%nodes(m)%total_weight
         if (den > 0.0_dp) then
            tree%nodes(idx)%prediction = (tree%nodes(l)%total_weight * tree%nodes(l)%prediction + &
                                          tree%nodes(r)%total_weight * tree%nodes(r)%prediction + &
                                          tree%nodes(m)%total_weight * tree%nodes(m)%prediction) / den
         end if
      end if
   end subroutine adjust_node

   subroutine adjust_tree(tree, min_obs, assignment, delta)
      type(gbm_tree), intent(inout) :: tree
      integer, intent(in) :: min_obs
      integer, intent(in) :: assignment(:)
      real(dp), intent(out) :: delta(size(assignment))
      integer :: i
      call adjust_node(tree, 1, min_obs)
      do i = 1, size(assignment)
         delta(i) = tree%nodes(assignment(i))%prediction
      end do
   end subroutine adjust_tree

   real(dp) function predict_tree(tree, row) result(pred)
      type(gbm_tree), intent(in) :: tree
      real(dp), intent(in) :: row(:)
      integer :: idx
      idx = 1
      do while (.not. tree%nodes(idx)%is_terminal)
         idx = choose_child(tree%nodes(idx), row)
      end do
      pred = tree%shrinkage * tree%nodes(idx)%prediction
   end function predict_tree

   recursive subroutine influence_node(tree, idx, influence)
      type(gbm_tree), intent(in) :: tree
      integer, intent(in) :: idx
      real(dp), intent(inout) :: influence(:)
      if (tree%nodes(idx)%is_terminal) return
      influence(tree%nodes(idx)%split_var) = influence(tree%nodes(idx)%split_var) + tree%nodes(idx)%improvement
      ! Preserve upstream behavior: influence traversal follows left/right but
      ! does not recurse through the missing branch.
      call influence_node(tree, tree%nodes(idx)%left, influence)
      call influence_node(tree, tree%nodes(idx)%right, influence)
   end subroutine influence_node

   subroutine tree_relative_influence(tree, influence)
      type(gbm_tree), intent(in) :: tree
      real(dp), intent(inout) :: influence(:)
      if (tree%n_nodes > 0) call influence_node(tree, 1, influence)
   end subroutine tree_relative_influence

end module gbm3_tree
