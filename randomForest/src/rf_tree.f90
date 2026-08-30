! SPDX-License-Identifier: GPL-2.0-or-later
module rf_tree_mod
   use r_kinds, only : dp
   use rf_rng, only : rf_rng_state, shuffle_int
   use rf_sort, only : sort_indices_by_values, sort_categories_by_values
   use rf_types, only : rf_tree, RF_TERMINAL, RF_INTERIOR
   implicit none
   private

   real(dp), parameter :: split_eps = 1.0e-12_dp

   public :: build_class_tree, predict_class_tree
   public :: build_reg_tree, predict_reg_tree, tree_var_used

contains

   subroutine initialize_tree(tree, maxnodes, maxcat, nvar)
      type(rf_tree), intent(out) :: tree
      integer, intent(in) :: maxnodes, maxcat, nvar

      tree%n_nodes = 0
      tree%maxcat = max(1, maxcat)
      allocate(tree%left(maxnodes), tree%right(maxnodes), tree%split_var(maxnodes))
      allocate(tree%status(maxnodes), tree%node_class(maxnodes))
      allocate(tree%split_value(maxnodes), tree%node_mean(maxnodes))
      allocate(tree%cat_left(tree%maxcat, maxnodes))
      allocate(tree%impurity_decrease(nvar))
      tree%left = 0
      tree%right = 0
      tree%split_var = 0
      tree%status = RF_TERMINAL
      tree%node_class = 0
      tree%split_value = 0.0_dp
      tree%node_mean = 0.0_dp
      tree%cat_left = .false.
      tree%impurity_decrease = 0.0_dp
   end subroutine initialize_tree

   subroutine build_class_tree(x, y, obs_index, obs_weight, ncat, nclass, mtry, nodesize, maxnodes, rng, tree)
      real(dp), intent(in) :: x(:,:)
      integer, intent(in) :: y(:), obs_index(:), ncat(:), nclass, mtry, nodesize, maxnodes
      real(dp), intent(in) :: obs_weight(:)
      type(rf_rng_state), intent(inout) :: rng
      type(rf_tree), intent(out) :: tree
      integer, allocatable :: node_start(:), node_count(:), work_index(:), temp_index(:)
      real(dp), allocatable :: work_weight(:), temp_weight(:), class_count(:)
      logical, allocatable :: best_mask(:)
      integer :: node, next_node, start, count_node, finish, i, j, nleft
      integer :: best_var, best_class, maxcat
      real(dp) :: best_cut, improve
      logical :: found

      maxcat = maxval(ncat)
      call initialize_tree(tree, maxnodes, maxcat, size(x, 2))
      allocate(node_start(maxnodes), node_count(maxnodes))
      allocate(work_index(size(obs_index)), temp_index(size(obs_index)))
      allocate(work_weight(size(obs_weight)), temp_weight(size(obs_weight)))
      allocate(class_count(nclass), best_mask(maxcat))
      node_start = 0
      node_count = 0
      work_index = obs_index
      work_weight = obs_weight
      node_start(1) = 1
      node_count(1) = size(obs_index)
      next_node = 1

      do node = 1, maxnodes
         if (node > next_node) exit
         start = node_start(node)
         count_node = node_count(node)
         if (count_node <= 0) cycle
         finish = start + count_node - 1

         class_count = 0.0_dp
         do i = start, finish
            class_count(y(work_index(i))) = class_count(y(work_index(i))) + work_weight(i)
         end do
         call choose_class(class_count, rng, best_class)
         tree%node_class(node) = best_class

         if (count_node <= nodesize .or. count(class_count > split_eps) <= 1) then
            tree%status(node) = RF_TERMINAL
            cycle
         end if
         if (next_node + 2 > maxnodes) then
            tree%status(node) = RF_TERMINAL
            cycle
         end if

         call best_class_split(x, y, work_index(start:finish), work_weight(start:finish), ncat, &
            nclass, mtry, rng, best_var, best_cut, best_mask, improve, found)
         if (.not. found .or. improve <= split_eps) then
            tree%status(node) = RF_TERMINAL
            cycle
         end if

         nleft = 0
         do i = start, finish
            j = work_index(i)
            if (go_left(x(j, best_var), ncat(best_var), best_cut, best_mask)) then
               nleft = nleft + 1
               temp_index(start + nleft - 1) = j
               temp_weight(start + nleft - 1) = work_weight(i)
            end if
         end do
         j = start + nleft
         do i = start, finish
            if (.not. go_left(x(work_index(i), best_var), ncat(best_var), best_cut, best_mask)) then
               temp_index(j) = work_index(i)
               temp_weight(j) = work_weight(i)
               j = j + 1
            end if
         end do
         if (nleft <= 0 .or. nleft >= count_node) then
            tree%status(node) = RF_TERMINAL
            cycle
         end if
         work_index(start:finish) = temp_index(start:finish)
         work_weight(start:finish) = temp_weight(start:finish)

         tree%status(node) = RF_INTERIOR
         tree%split_var(node) = best_var
         tree%split_value(node) = best_cut
         if (ncat(best_var) > 1) tree%cat_left(1:ncat(best_var), node) = best_mask(1:ncat(best_var))
         tree%impurity_decrease(best_var) = tree%impurity_decrease(best_var) + max(improve, 0.0_dp)

         tree%left(node) = next_node + 1
         tree%right(node) = next_node + 2
         node_start(next_node + 1) = start
         node_count(next_node + 1) = nleft
         node_start(next_node + 2) = start + nleft
         node_count(next_node + 2) = count_node - nleft
         next_node = next_node + 2
      end do
      tree%n_nodes = next_node
   end subroutine build_class_tree

   subroutine best_class_split(x, y, indices, weights, ncat, nclass, mtry, rng, &
      best_var, best_cut, best_mask, best_improve, found)
      real(dp), intent(in) :: x(:,:), weights(:)
      integer, intent(in) :: y(:), indices(:), ncat(:), nclass, mtry
      type(rf_rng_state), intent(inout) :: rng
      integer, intent(out) :: best_var
      real(dp), intent(out) :: best_cut, best_improve
      logical, intent(out) :: best_mask(:), found
      integer, allocatable :: vars(:)
      logical, allocatable :: candidate_mask(:)
      real(dp), allocatable :: parent_count(:)
      integer :: i, v, ntry
      real(dp) :: parent_score, crit, cut
      logical :: ok

      allocate(vars(size(x, 2)), candidate_mask(size(best_mask)), parent_count(nclass))
      vars = [(i, i = 1, size(vars))]
      call shuffle_int(rng, vars)
      parent_count = 0.0_dp
      do i = 1, size(indices)
         parent_count(y(indices(i))) = parent_count(y(indices(i))) + weights(i)
      end do
      parent_score = sum(parent_count * parent_count) / max(sum(parent_count), split_eps)

      best_var = 0
      best_cut = 0.0_dp
      best_mask = .false.
      best_improve = -huge(1.0_dp)
      found = .false.
      ntry = min(max(1, mtry), size(vars))

      do i = 1, ntry
         v = vars(i)
         candidate_mask = .false.
         if (ncat(v) <= 1) then
            call class_numeric_split(x(:, v), y, indices, weights, nclass, rng, crit, cut, ok)
         else
            call class_categorical_split(x(:, v), y, indices, weights, nclass, ncat(v), rng, &
               crit, candidate_mask, ok)
            cut = 0.0_dp
         end if
         if (.not. ok) cycle
         if (.not. found .or. crit - parent_score > best_improve + tie_tolerance(crit, parent_score)) then
            found = .true.
            best_var = v
            best_cut = cut
            best_mask = candidate_mask
            best_improve = crit - parent_score
         else if (is_tie(crit - parent_score, best_improve)) then
            if (rng%uniform() < 0.5_dp) then
               best_var = v
               best_cut = cut
               best_mask = candidate_mask
               best_improve = crit - parent_score
            end if
         end if
      end do
   end subroutine best_class_split

   subroutine class_numeric_split(values, y, indices, weights, nclass, rng, best_crit, best_cut, found)
      real(dp), intent(in) :: values(:), weights(:)
      integer, intent(in) :: y(:), indices(:), nclass
      type(rf_rng_state), intent(inout) :: rng
      real(dp), intent(out) :: best_crit, best_cut
      logical, intent(out) :: found
      integer, allocatable :: order(:)
      real(dp), allocatable :: left_count(:), right_count(:), weight_by_obs(:)
      integer :: i, obs, next_obs, ntie
      real(dp) :: left_den, right_den, crit, midpoint

      allocate(order(size(indices)), left_count(nclass), right_count(nclass), weight_by_obs(size(values)))
      order = indices
      weight_by_obs = 0.0_dp
      do i = 1, size(indices)
         weight_by_obs(indices(i)) = weights(i)
      end do
      call sort_indices_by_values(values, order)
      left_count = 0.0_dp
      right_count = 0.0_dp
      do i = 1, size(indices)
         right_count(y(indices(i))) = right_count(y(indices(i))) + weights(i)
      end do
      left_den = 0.0_dp
      right_den = sum(weights)
      best_crit = -huge(1.0_dp)
      best_cut = 0.0_dp
      found = .false.
      ntie = 1

      do i = 1, size(order) - 1
         obs = order(i)
         next_obs = order(i + 1)
         left_count(y(obs)) = left_count(y(obs)) + weight_by_obs(obs)
         right_count(y(obs)) = right_count(y(obs)) - weight_by_obs(obs)
         left_den = left_den + weight_by_obs(obs)
         right_den = right_den - weight_by_obs(obs)
         if (values(obs) >= values(next_obs) - split_eps) cycle
         if (min(left_den, right_den) <= 1.0e-8_dp) cycle
         crit = sum(left_count * left_count) / left_den + sum(right_count * right_count) / right_den
         midpoint = 0.5_dp * (values(obs) + values(next_obs))
         if (.not. found .or. crit > best_crit + tie_tolerance(crit, best_crit)) then
            found = .true.
            best_crit = crit
            best_cut = midpoint
            ntie = 1
         else if (is_tie(crit, best_crit)) then
            ntie = ntie + 1
            if (rng%uniform() < 1.0_dp / real(ntie, dp)) best_cut = midpoint
         end if
      end do
   end subroutine class_numeric_split

   subroutine class_categorical_split(values, y, indices, weights, nclass, ncat, rng, best_crit, best_mask, found)
      real(dp), intent(in) :: values(:), weights(:)
      integer, intent(in) :: y(:), indices(:), nclass, ncat
      type(rf_rng_state), intent(inout) :: rng
      real(dp), intent(out) :: best_crit
      logical, intent(out) :: best_mask(:), found
      real(dp), allocatable :: table(:,:), class_total(:), cat_total(:), left(:), proportions(:)
      integer, allocatable :: categories(:)
      logical, allocatable :: mask(:)
      integer :: i, j, obs, cat_value, nsplit, pattern, ntie
      real(dp) :: total_weight, left_den, right_den, crit

      allocate(table(nclass, ncat), class_total(nclass), cat_total(ncat), left(nclass))
      allocate(proportions(ncat), categories(ncat), mask(ncat))
      table = 0.0_dp
      class_total = 0.0_dp
      cat_total = 0.0_dp
      do i = 1, size(indices)
         obs = indices(i)
         cat_value = nint(values(obs))
         if (cat_value < 1 .or. cat_value > ncat) cycle
         table(y(obs), cat_value) = table(y(obs), cat_value) + weights(i)
         class_total(y(obs)) = class_total(y(obs)) + weights(i)
         cat_total(cat_value) = cat_total(cat_value) + weights(i)
      end do
      total_weight = sum(class_total)
      found = .false.
      best_crit = -huge(1.0_dp)
      best_mask = .false.
      ntie = 1
      if (count(cat_total > split_eps) <= 1) return

      if (nclass == 2 .and. ncat > 10) then
         proportions = 0.0_dp
         do i = 1, ncat
            if (cat_total(i) > split_eps) proportions(i) = table(1, i) / cat_total(i)
         end do
         call sort_categories_by_values(proportions, categories, ncat)
         left = 0.0_dp
         left_den = 0.0_dp
         mask = .false.
         do i = 1, ncat - 1
            j = categories(i)
            left = left + table(:, j)
            left_den = left_den + cat_total(j)
            mask(j) = .true.
            right_den = total_weight - left_den
            if (proportions(j) >= proportions(categories(i + 1)) - split_eps) cycle
            if (min(left_den, right_den) <= 1.0e-8_dp) cycle
            crit = sum(left * left) / left_den + sum((class_total - left) ** 2) / right_den
            call update_cat_candidate(crit, mask, rng, best_crit, best_mask, found, ntie)
         end do
         return
      end if

      if (ncat <= 10) then
         nsplit = 2 ** (ncat - 1) - 1
         do pattern = 1, nsplit
            do j = 1, ncat
               mask(j) = btest(pattern, j - 1)
            end do
            if (.not. valid_cat_mask(mask, cat_total)) cycle
            call evaluate_cat_mask(mask, table, class_total, cat_total, total_weight, crit)
            call update_cat_candidate(crit, mask, rng, best_crit, best_mask, found, ntie)
         end do
      else
         nsplit = 512
         do pattern = 1, nsplit
            do j = 1, ncat
               mask(j) = rng%uniform() > 0.5_dp
            end do
            if (.not. valid_cat_mask(mask, cat_total)) cycle
            call evaluate_cat_mask(mask, table, class_total, cat_total, total_weight, crit)
            call update_cat_candidate(crit, mask, rng, best_crit, best_mask, found, ntie)
         end do
      end if
   end subroutine class_categorical_split

   subroutine evaluate_cat_mask(mask, table, class_total, cat_total, total_weight, crit, found_mask)
      logical, intent(in) :: mask(:)
      real(dp), intent(in) :: table(:,:), class_total(:), cat_total(:), total_weight
      real(dp), intent(out) :: crit
      logical, intent(out), optional :: found_mask
      real(dp), allocatable :: left(:)
      real(dp) :: left_den, right_den
      integer :: j

      allocate(left(size(class_total)))
      left = 0.0_dp
      left_den = 0.0_dp
      do j = 1, size(mask)
         if (mask(j)) then
            left = left + table(:, j)
            left_den = left_den + cat_total(j)
         end if
      end do
      right_den = total_weight - left_den
      if (present(found_mask)) found_mask = min(left_den, right_den) > 1.0e-8_dp
      if (min(left_den, right_den) <= 1.0e-8_dp) then
         crit = -huge(1.0_dp)
         return
      end if
      crit = sum(left * left) / left_den + sum((class_total - left) ** 2) / right_den
   end subroutine evaluate_cat_mask

   logical function valid_cat_mask(mask, cat_total)
      logical, intent(in) :: mask(:)
      real(dp), intent(in) :: cat_total(:)
      real(dp) :: left_den, total

      left_den = sum(cat_total, mask=mask)
      total = sum(cat_total)
      valid_cat_mask = left_den > 1.0e-8_dp .and. total - left_den > 1.0e-8_dp
   end function valid_cat_mask

   subroutine update_cat_candidate(crit, mask, rng, best_crit, best_mask, found, ntie)
      real(dp), intent(in) :: crit
      logical, intent(in) :: mask(:)
      type(rf_rng_state), intent(inout) :: rng
      real(dp), intent(inout) :: best_crit
      logical, intent(inout) :: best_mask(:), found
      integer, intent(inout) :: ntie

      if (.not. found .or. crit > best_crit + tie_tolerance(crit, best_crit)) then
         found = .true.
         best_crit = crit
         best_mask(1:size(mask)) = mask
         ntie = 1
      else if (is_tie(crit, best_crit)) then
         ntie = ntie + 1
         if (rng%uniform() < 1.0_dp / real(ntie, dp)) best_mask(1:size(mask)) = mask
      end if
   end subroutine update_cat_candidate

   subroutine choose_class(counts, rng, klass)
      real(dp), intent(in) :: counts(:)
      type(rf_rng_state), intent(inout) :: rng
      integer, intent(out) :: klass
      integer :: j, ntie
      real(dp) :: best

      klass = 1
      best = counts(1)
      ntie = 1
      do j = 2, size(counts)
         if (counts(j) > best + tie_tolerance(counts(j), best)) then
            best = counts(j)
            klass = j
            ntie = 1
         else if (is_tie(counts(j), best)) then
            ntie = ntie + 1
            if (rng%uniform() < 1.0_dp / real(ntie, dp)) klass = j
         end if
      end do
   end subroutine choose_class

   subroutine predict_class_tree(tree, x, ncat, prediction, terminal_node)
      type(rf_tree), intent(in) :: tree
      real(dp), intent(in) :: x(:,:)
      integer, intent(in) :: ncat(:)
      integer, intent(out) :: prediction(size(x, 1))
      integer, intent(out), optional :: terminal_node(size(x, 1))
      integer :: i, node, v

      do i = 1, size(x, 1)
         node = 1
         do while (node <= tree%n_nodes .and. tree%status(node) == RF_INTERIOR)
            v = tree%split_var(node)
            if (ncat(v) <= 1) then
               if (x(i, v) <= tree%split_value(node)) then
                  node = tree%left(node)
               else
                  node = tree%right(node)
               end if
            else
               if (category_goes_left(x(i, v), ncat(v), tree%cat_left(:, node))) then
                  node = tree%left(node)
               else
                  node = tree%right(node)
               end if
            end if
         end do
         if (node < 1 .or. node > tree%n_nodes) node = 1
         prediction(i) = tree%node_class(node)
         if (present(terminal_node)) terminal_node(i) = node
      end do
   end subroutine predict_class_tree

   subroutine build_reg_tree(x, y, obs_index, obs_weight, ncat, mtry, nodesize, maxnodes, rng, tree)
      real(dp), intent(in) :: x(:,:), y(:), obs_weight(:)
      integer, intent(in) :: obs_index(:), ncat(:), mtry, nodesize, maxnodes
      type(rf_rng_state), intent(inout) :: rng
      type(rf_tree), intent(out) :: tree
      integer, allocatable :: node_start(:), node_count(:), work_index(:), temp_index(:)
      real(dp), allocatable :: work_weight(:), temp_weight(:)
      logical, allocatable :: best_mask(:)
      integer :: node, next_node, start, finish, count_node, i, j, nleft, best_var, maxcat
      real(dp) :: best_cut, improve, sw, sy
      logical :: found

      maxcat = maxval(ncat)
      call initialize_tree(tree, maxnodes, maxcat, size(x, 2))
      allocate(node_start(maxnodes), node_count(maxnodes))
      allocate(work_index(size(obs_index)), temp_index(size(obs_index)))
      allocate(work_weight(size(obs_weight)), temp_weight(size(obs_weight)))
      allocate(best_mask(maxcat))
      node_start = 0
      node_count = 0
      work_index = obs_index
      work_weight = obs_weight
      node_start(1) = 1
      node_count(1) = size(obs_index)
      next_node = 1

      do node = 1, maxnodes
         if (node > next_node) exit
         start = node_start(node)
         count_node = node_count(node)
         if (count_node <= 0) cycle
         finish = start + count_node - 1
         sw = sum(work_weight(start:finish))
         sy = 0.0_dp
         do i = start, finish
            sy = sy + work_weight(i) * y(work_index(i))
         end do
         tree%node_mean(node) = sy / max(sw, split_eps)

         if (sw <= real(nodesize, dp) + split_eps .or. next_node + 2 > maxnodes) then
            tree%status(node) = RF_TERMINAL
            cycle
         end if

         call best_reg_split(x, y, work_index(start:finish), work_weight(start:finish), ncat, &
            mtry, rng, best_var, best_cut, best_mask, improve, found)
         if (.not. found .or. improve <= split_eps) then
            tree%status(node) = RF_TERMINAL
            cycle
         end if

         nleft = 0
         do i = start, finish
            j = work_index(i)
            if (go_left(x(j, best_var), ncat(best_var), best_cut, best_mask)) then
               nleft = nleft + 1
               temp_index(start + nleft - 1) = j
               temp_weight(start + nleft - 1) = work_weight(i)
            end if
         end do
         j = start + nleft
         do i = start, finish
            if (.not. go_left(x(work_index(i), best_var), ncat(best_var), best_cut, best_mask)) then
               temp_index(j) = work_index(i)
               temp_weight(j) = work_weight(i)
               j = j + 1
            end if
         end do
         if (nleft <= 0 .or. nleft >= count_node) then
            tree%status(node) = RF_TERMINAL
            cycle
         end if
         work_index(start:finish) = temp_index(start:finish)
         work_weight(start:finish) = temp_weight(start:finish)

         tree%status(node) = RF_INTERIOR
         tree%split_var(node) = best_var
         tree%split_value(node) = best_cut
         if (ncat(best_var) > 1) tree%cat_left(1:ncat(best_var), node) = best_mask(1:ncat(best_var))
         tree%impurity_decrease(best_var) = tree%impurity_decrease(best_var) + improve
         tree%left(node) = next_node + 1
         tree%right(node) = next_node + 2
         node_start(next_node + 1) = start
         node_count(next_node + 1) = nleft
         node_start(next_node + 2) = start + nleft
         node_count(next_node + 2) = count_node - nleft
         next_node = next_node + 2
      end do
      tree%n_nodes = next_node
   end subroutine build_reg_tree

   subroutine best_reg_split(x, y, indices, weights, ncat, mtry, rng, best_var, best_cut, best_mask, &
      best_improve, found)
      real(dp), intent(in) :: x(:,:), y(:), weights(:)
      integer, intent(in) :: indices(:), ncat(:), mtry
      type(rf_rng_state), intent(inout) :: rng
      integer, intent(out) :: best_var
      real(dp), intent(out) :: best_cut, best_improve
      logical, intent(out) :: best_mask(:), found
      integer, allocatable :: vars(:)
      logical, allocatable :: candidate_mask(:)
      integer :: i, v, ntry
      real(dp) :: crit, cut
      logical :: ok

      allocate(vars(size(x, 2)), candidate_mask(size(best_mask)))
      vars = [(i, i = 1, size(vars))]
      call shuffle_int(rng, vars)
      best_var = 0
      best_cut = 0.0_dp
      best_mask = .false.
      best_improve = -huge(1.0_dp)
      found = .false.
      ntry = min(max(1, mtry), size(vars))

      do i = 1, ntry
         v = vars(i)
         candidate_mask = .false.
         if (ncat(v) <= 1) then
            call reg_numeric_split(x(:, v), y, indices, weights, rng, crit, cut, ok)
         else
            call reg_categorical_split(x(:, v), y, indices, weights, ncat(v), rng, crit, candidate_mask, ok)
            cut = 0.0_dp
         end if
         if (.not. ok) cycle
         if (.not. found .or. crit > best_improve + tie_tolerance(crit, best_improve)) then
            found = .true.
            best_var = v
            best_cut = cut
            best_mask = candidate_mask
            best_improve = crit
         else if (is_tie(crit, best_improve)) then
            if (rng%uniform() < 0.5_dp) then
               best_var = v
               best_cut = cut
               best_mask = candidate_mask
               best_improve = crit
            end if
         end if
      end do
   end subroutine best_reg_split

   subroutine reg_numeric_split(values, y, indices, weights, rng, best_crit, best_cut, found)
      real(dp), intent(in) :: values(:), y(:), weights(:)
      integer, intent(in) :: indices(:)
      type(rf_rng_state), intent(inout) :: rng
      real(dp), intent(out) :: best_crit, best_cut
      logical, intent(out) :: found
      integer, allocatable :: order(:)
      real(dp), allocatable :: weight_by_obs(:)
      integer :: i, obs, next_obs, ntie
      real(dp) :: total_w, total_y, left_w, left_y, right_w, right_y, parent, crit, midpoint

      allocate(order(size(indices)), weight_by_obs(size(values)))
      order = indices
      weight_by_obs = 0.0_dp
      do i = 1, size(indices)
         weight_by_obs(indices(i)) = weights(i)
      end do
      call sort_indices_by_values(values, order)
      total_w = sum(weights)
      total_y = 0.0_dp
      do i = 1, size(indices)
         total_y = total_y + weights(i) * y(indices(i))
      end do
      parent = total_y * total_y / max(total_w, split_eps)
      left_w = 0.0_dp
      left_y = 0.0_dp
      best_crit = -huge(1.0_dp)
      best_cut = 0.0_dp
      found = .false.
      ntie = 1

      do i = 1, size(order) - 1
         obs = order(i)
         next_obs = order(i + 1)
         left_w = left_w + weight_by_obs(obs)
         left_y = left_y + weight_by_obs(obs) * y(obs)
         right_w = total_w - left_w
         right_y = total_y - left_y
         if (values(obs) >= values(next_obs) - split_eps) cycle
         if (min(left_w, right_w) <= 1.0e-8_dp) cycle
         crit = left_y * left_y / left_w + right_y * right_y / right_w - parent
         midpoint = 0.5_dp * (values(obs) + values(next_obs))
         if (.not. found .or. crit > best_crit + tie_tolerance(crit, best_crit)) then
            found = .true.
            best_crit = crit
            best_cut = midpoint
            ntie = 1
         else if (is_tie(crit, best_crit)) then
            ntie = ntie + 1
            if (rng%uniform() < 1.0_dp / real(ntie, dp)) best_cut = midpoint
         end if
      end do
   end subroutine reg_numeric_split

   subroutine reg_categorical_split(values, y, indices, weights, ncat, rng, best_crit, best_mask, found)
      real(dp), intent(in) :: values(:), y(:), weights(:)
      integer, intent(in) :: indices(:), ncat
      type(rf_rng_state), intent(inout) :: rng
      real(dp), intent(out) :: best_crit
      logical, intent(out) :: best_mask(:), found
      real(dp), allocatable :: cat_w(:), cat_y(:), cat_mean(:)
      integer, allocatable :: categories(:)
      logical, allocatable :: mask(:)
      integer :: i, obs, c, ntie
      real(dp) :: total_w, total_y, left_w, left_y, right_w, right_y, parent, crit

      allocate(cat_w(ncat), cat_y(ncat), cat_mean(ncat), categories(ncat), mask(ncat))
      cat_w = 0.0_dp
      cat_y = 0.0_dp
      do i = 1, size(indices)
         obs = indices(i)
         c = nint(values(obs))
         if (c < 1 .or. c > ncat) cycle
         cat_w(c) = cat_w(c) + weights(i)
         cat_y(c) = cat_y(c) + weights(i) * y(obs)
      end do
      cat_mean = 0.0_dp
      do c = 1, ncat
         if (cat_w(c) > split_eps) cat_mean(c) = cat_y(c) / cat_w(c)
      end do
      call sort_categories_by_values(cat_mean, categories, ncat)
      total_w = sum(cat_w)
      total_y = sum(cat_y)
      parent = total_y * total_y / max(total_w, split_eps)
      left_w = 0.0_dp
      left_y = 0.0_dp
      mask = .false.
      best_mask = .false.
      best_crit = -huge(1.0_dp)
      found = .false.
      ntie = 1
      do i = 1, ncat - 1
         c = categories(i)
         left_w = left_w + cat_w(c)
         left_y = left_y + cat_y(c)
         mask(c) = .true.
         right_w = total_w - left_w
         right_y = total_y - left_y
         if (cat_mean(c) >= cat_mean(categories(i + 1)) - split_eps) cycle
         if (min(left_w, right_w) <= 1.0e-8_dp) cycle
         crit = left_y * left_y / left_w + right_y * right_y / right_w - parent
         if (.not. found .or. crit > best_crit + tie_tolerance(crit, best_crit)) then
            found = .true.
            best_crit = crit
            best_mask(1:ncat) = mask
            ntie = 1
         else if (is_tie(crit, best_crit)) then
            ntie = ntie + 1
            if (rng%uniform() < 1.0_dp / real(ntie, dp)) best_mask(1:ncat) = mask
         end if
      end do
   end subroutine reg_categorical_split

   subroutine predict_reg_tree(tree, x, ncat, prediction, terminal_node)
      type(rf_tree), intent(in) :: tree
      real(dp), intent(in) :: x(:,:)
      integer, intent(in) :: ncat(:)
      real(dp), intent(out) :: prediction(size(x, 1))
      integer, intent(out), optional :: terminal_node(size(x, 1))
      integer :: i, node, v

      do i = 1, size(x, 1)
         node = 1
         do while (node <= tree%n_nodes .and. tree%status(node) == RF_INTERIOR)
            v = tree%split_var(node)
            if (ncat(v) <= 1) then
               if (x(i, v) <= tree%split_value(node)) then
                  node = tree%left(node)
               else
                  node = tree%right(node)
               end if
            else
               if (category_goes_left(x(i, v), ncat(v), tree%cat_left(:, node))) then
                  node = tree%left(node)
               else
                  node = tree%right(node)
               end if
            end if
         end do
         if (node < 1 .or. node > tree%n_nodes) node = 1
         prediction(i) = tree%node_mean(node)
         if (present(terminal_node)) terminal_node(i) = node
      end do
   end subroutine predict_reg_tree

   logical function go_left(value, ncat, cut, mask)
      real(dp), intent(in) :: value, cut
      integer, intent(in) :: ncat
      logical, intent(in) :: mask(:)

      if (ncat <= 1) then
         go_left = value <= cut
      else
         go_left = category_goes_left(value, ncat, mask)
      end if
   end function go_left

   logical function category_goes_left(value, ncat, mask)
      real(dp), intent(in) :: value
      integer, intent(in) :: ncat
      logical, intent(in) :: mask(:)
      integer :: c

      c = nint(value)
      if (c < 1 .or. c > ncat) then
         category_goes_left = .false.
      else
         category_goes_left = mask(c)
      end if
   end function category_goes_left

   subroutine tree_var_used(tree, used)
      type(rf_tree), intent(in) :: tree
      logical, intent(out) :: used(:)
      integer :: node, v

      used = .false.
      do node = 1, tree%n_nodes
         if (tree%status(node) == RF_INTERIOR) then
            v = tree%split_var(node)
            if (v >= 1 .and. v <= size(used)) used(v) = .true.
         end if
      end do
   end subroutine tree_var_used

   pure logical function is_tie(a, b)
      real(dp), intent(in) :: a, b
      is_tie = abs(a - b) <= tie_tolerance(a, b)
   end function is_tie

   pure real(dp) function tie_tolerance(a, b)
      real(dp), intent(in) :: a, b
      tie_tolerance = 64.0_dp * epsilon(1.0_dp) * max(1.0_dp, abs(a), abs(b))
   end function tie_tolerance

end module rf_tree_mod
