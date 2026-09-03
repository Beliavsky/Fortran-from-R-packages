! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
!
! Triangle reconstruction follows ape src/triangMtd.c and src/triangMtds.c
! (Copyright 2011-2012 Andrei-Alin Popescu).
module ape_triangle_reconstruction
   use ieee_arithmetic, only : ieee_is_nan
   use r_kinds, only : dp
   use ape_types, only : phylo_tree, make_phylo_tree
   implicit none
   private

   public :: triang_mtd
   public :: triang_mtds

contains

   pure subroutine triang_mtd(distance, tree, info)
      !! Reconstructs a tree from complete distances with ape's triangles method.
      real(dp), intent(in) :: distance(:, :) !! Symmetric complete taxon-distance matrix; diagonal entries are ignored.
      type(phylo_tree), intent(out) :: tree !! Reconstructed ape-style rooted encoding of the unrooted inferred tree.
      integer, intent(out) :: info !! Status code: zero on success; nonzero indicates invalid distances or attachment geometry.
      integer, allocatable :: edge(:, :)
      integer, allocatable :: edge_start(:)
      integer, allocatable :: edge_end(:)
      integer, allocatable :: path(:)
      real(dp), allocatable :: edge_length(:)
      real(dp), allocatable :: leaf_distance(:)
      logical, allocatable :: present(:)
      integer :: attachment_edge
      integer :: aux
      integer :: child
      integer :: edge_count
      integer :: i
      integer :: j
      integer :: n
      integer :: new_node
      integer :: path_index
      integer :: present_count
      integer :: x
      integer :: x3
      integer :: y
      integer :: z
      real(dp) :: attach_from_x
      real(dp) :: minimum_distance
      real(dp) :: new_distance
      real(dp) :: previous_sum
      real(dp) :: segment_length
      real(dp) :: path_sum
      logical :: forward

      info = 0
      n = size(distance, 1)
      if (n < 3 .or. size(distance, 2) /= n) then
         info = 1
         return
      end if
      if (.not. valid_complete_distance_matrix(distance)) then
         info = 2
         return
      end if

      allocate(edge(2 * n - 3, 2), edge_length(2 * n - 3))
      allocate(edge_start(n), edge_end(n), leaf_distance(n), present(n))
      edge = 0
      edge_length = 0.0_dp
      edge_start = 0
      edge_end = 0
      leaf_distance = huge(1.0_dp)
      present = .false.

      call minimum_positive_pair(distance, x, y, info)
      if (info /= 0) return
      present(x) = .true.
      present(y) = .true.
      do i = 1, n
         if (present(i)) cycle
         edge_start(i) = x
         edge_end(i) = y
         leaf_distance(i) = 0.5_dp * (distance(i, x) + distance(i, y) - distance(x, y))
      end do

      minimum_distance = huge(1.0_dp)
      x3 = 0
      do i = 1, n
         if (present(i)) cycle
         if (leaf_distance(i) < minimum_distance) then
            minimum_distance = leaf_distance(i)
            x3 = i
         end if
      end do
      if (x3 == 0) then
         info = 3
         return
      end if

      new_node = n + 1
      edge_count = 3
      edge(1, :) = [new_node, x]
      edge_length(1) = 0.5_dp * (distance(x3, x) + distance(y, x) - distance(y, x3))
      edge(2, :) = [new_node, y]
      edge_length(2) = 0.5_dp * (distance(y, x3) + distance(y, x) - distance(x, x3))
      edge(3, :) = [new_node, x3]
      edge_length(3) = 0.5_dp * (distance(x3, x) + distance(y, x3) - distance(y, x))
      present(x3) = .true.
      present_count = 3

      do z = 1, n
         if (present(z)) cycle
         do i = 1, n
            if (i == x3 .or. .not. present(i)) cycle
            new_distance = 0.5_dp * (distance(i, z) + distance(z, x3) - distance(i, x3))
            if (new_distance < leaf_distance(z)) then
               leaf_distance(z) = new_distance
               edge_start(z) = i
               edge_end(z) = x3
            end if
         end do
      end do

      do while (present_count < n)
         minimum_distance = huge(1.0_dp)
         z = 0
         do i = 1, n
            if (present(i)) cycle
            if (leaf_distance(i) < minimum_distance) then
               minimum_distance = leaf_distance(i)
               z = i
            end if
         end do
         if (z == 0) then
            info = 4
            return
         end if
         x = edge_start(z)
         y = edge_end(z)
         call partial_tree_path(x, y, n + n - 2, edge, edge_count, path, info)
         if (info /= 0) return

         attach_from_x = 0.5_dp * (distance(y, x) + distance(z, x) - distance(z, y))
         path_sum = 0.0_dp
         previous_sum = 0.0_dp
         attachment_edge = 0
         forward = .false.
         do path_index = 1, size(path) - 1
            aux = path(path_index)
            child = path(path_index + 1)
            call partial_edge_lookup(aux, child, edge, edge_length, edge_count, attachment_edge, &
               segment_length, forward)
            if (attachment_edge == 0) then
               info = 5
               return
            end if
            previous_sum = path_sum
            path_sum = path_sum + segment_length
            if (path_sum >= attach_from_x) exit
         end do
         if (attachment_edge == 0 .or. path_sum < attach_from_x) then
            info = 5
            return
         end if

         new_node = new_node + 1
         child = edge(attachment_edge, 2)
         edge(attachment_edge, 2) = new_node
         if (forward) then
            edge_length(attachment_edge) = attach_from_x - previous_sum
         else
            edge_length(attachment_edge) = path_sum - attach_from_x
         end if
         edge_count = edge_count + 1
         edge(edge_count, :) = [new_node, child]
         if (forward) then
            edge_length(edge_count) = path_sum - attach_from_x
         else
            edge_length(edge_count) = attach_from_x - previous_sum
         end if
         edge_count = edge_count + 1
         edge(edge_count, :) = [new_node, z]
         edge_length(edge_count) = minimum_distance

         present(z) = .true.
         present_count = present_count + 1
         do j = 1, n
            if (present(j)) cycle
            do i = 1, n
               if (i == z) cycle
               if (i /= x .and. i /= y) cycle
               new_distance = 0.5_dp * (distance(i, j) + distance(z, j) - distance(i, z))
               if (new_distance < leaf_distance(j)) then
                  leaf_distance(j) = new_distance
                  edge_start(j) = i
                  edge_end(j) = z
               end if
            end do
         end do
      end do

      tree = make_phylo_tree(n, edge(1:edge_count, :), edge_length(1:edge_count))
   end subroutine triang_mtd

   pure subroutine triang_mtds(distance, tree, info)
      !! Reconstructs a tree from incomplete distances with ape's sparse triangles method.
      real(dp), intent(in) :: distance(:, :) !! Symmetric matrix; NaN or negative off-diagonal entries denote missing distances.
      type(phylo_tree), intent(out) :: tree !! Tree from greedy complete-subset initialization and attachment.
      integer, intent(out) :: info !! Status code: zero on success; nonzero indicates invalid or insufficient distance information.
      real(dp), allocatable :: d(:, :)
      real(dp), allocatable :: sub_distance(:, :)
      integer, allocatable :: complete_order(:)
      integer, allocatable :: known_count(:)
      integer, allocatable :: old_tip(:)
      integer, allocatable :: edge(:, :)
      real(dp), allocatable :: edge_length(:)
      integer, allocatable :: path(:)
      logical, allocatable :: excluded(:)
      type(phylo_tree) :: seed
      integer :: attachment_edge
      integer :: aux
      integer :: child
      integer :: edge_count
      integer :: i
      integer :: ii
      integer :: j
      integer :: k
      integer :: max_count
      integer :: max_position
      integer :: n
      integer :: new_node
      integer :: path_index
      integer :: seed_start
      integer :: x
      integer :: y
      integer :: z
      real(dp) :: attach_from_x
      real(dp) :: minimum_distance
      real(dp) :: path_sum
      real(dp) :: previous_sum
      real(dp) :: segment_length
      real(dp) :: triangle_distance
      logical :: forward

      info = 0
      n = size(distance, 1)
      if (n < 3 .or. size(distance, 2) /= n) then
         info = 1
         return
      end if
      if (.not. valid_incomplete_distance_matrix_local(distance)) then
         info = 2
         return
      end if
      allocate(d(n, n))
      d = distance
      do i = 1, n
         d(i, i) = 0.0_dp
         do j = i + 1, n
            if (missing_distance(d(i, j))) then
               d(i, j) = -1.0_dp
               d(j, i) = -1.0_dp
            end if
         end do
      end do
      if (all_distances_known(d)) then
         call triang_mtd(d, tree, info)
         return
      end if

      allocate(known_count(n), complete_order(n), excluded(n))
      do i = 1, n
         known_count(i) = 1
         complete_order(i) = i
         do j = 1, n
            if (i == j) cycle
            if (.not. missing_distance(d(i, j))) known_count(i) = known_count(i) + 1
         end do
      end do
      call sort_by_known_count(known_count, complete_order)
      excluded = .false.
      seed_start = 1
      do while (seed_start <= n .and. known_count(seed_start) == n)
         seed_start = seed_start + 1
      end do
      do i = seed_start, n - 1
         if (excluded(complete_order(i))) cycle
         do j = i + 1, n
            if (excluded(complete_order(j))) cycle
            if (missing_distance(d(complete_order(i), complete_order(j)))) excluded(complete_order(j)) = .true.
         end do
      end do

      k = count(.not. excluded)
      if (k < 3) then
         info = 3
         return
      end if
      allocate(old_tip(k), sub_distance(k, k))
      ii = 0
      do i = 1, n
         if (excluded(i)) cycle
         ii = ii + 1
         old_tip(ii) = i
      end do
      do i = 1, k
         do j = 1, k
            sub_distance(i, j) = d(old_tip(i), old_tip(j))
         end do
      end do
      call triang_mtd(sub_distance, seed, info)
      if (info /= 0) return

      allocate(edge(2 * n - 3, 2), edge_length(2 * n - 3))
      edge = 0
      edge_length = 0.0_dp
      edge_count = seed%nedge()
      do i = 1, edge_count
         do j = 1, 2
            if (seed%edge(i, j) > k) then
               edge(i, j) = seed%edge(i, j) + n - k
            else
               edge(i, j) = old_tip(seed%edge(i, j))
            end if
         end do
         edge_length(i) = seed%edge_length(i)
      end do

      known_count = 0
      do i = 1, n
         if (.not. excluded(i)) cycle
         do j = 1, n
            if (excluded(j)) cycle
            if (.not. missing_distance(d(i, j))) known_count(i) = known_count(i) + 1
         end do
      end do
      new_node = n + k - 2

      do while (k < n)
         max_count = -1
         max_position = 0
         do i = 1, n
            if (.not. excluded(i)) cycle
            if (known_count(i) > max_count) then
               max_count = known_count(i)
               max_position = i
            end if
         end do
         if (max_position == 0) then
            info = 4
            return
         end if
         z = max_position
         excluded(z) = .false.
         do i = 1, n
            if (.not. excluded(i)) cycle
            if (.not. missing_distance(d(i, z))) known_count(i) = known_count(i) + 1
         end do

         minimum_distance = huge(1.0_dp)
         x = 0
         y = 0
         do i = 1, n - 1
            if (excluded(i) .or. i == z .or. missing_distance(d(i, z))) cycle
            do j = i + 1, n
               if (excluded(j) .or. j == z .or. missing_distance(d(j, z))) cycle
               triangle_distance = 0.5_dp * (d(i, z) + d(j, z) - d(i, j))
               if (triangle_distance < minimum_distance) then
                  minimum_distance = triangle_distance
                  x = i
                  y = j
               end if
            end do
         end do
         if (x == 0 .or. y == 0) then
            info = 5
            return
         end if

         call partial_tree_path(x, y, 2 * n - 2, edge, edge_count, path, info)
         if (info /= 0) return
         attach_from_x = 0.5_dp * (d(y, x) + d(z, x) - d(z, y))
         path_sum = 0.0_dp
         previous_sum = 0.0_dp
         attachment_edge = 0
         forward = .false.
         do path_index = 1, size(path) - 1
            aux = path(path_index)
            child = path(path_index + 1)
            call partial_edge_lookup(aux, child, edge, edge_length, edge_count, attachment_edge, &
               segment_length, forward)
            if (attachment_edge == 0) then
               info = 6
               return
            end if
            previous_sum = path_sum
            path_sum = path_sum + segment_length
            if (path_sum >= attach_from_x) exit
         end do
         if (attachment_edge == 0 .or. path_sum < attach_from_x) then
            info = 6
            return
         end if

         new_node = new_node + 1
         child = edge(attachment_edge, 2)
         edge(attachment_edge, 2) = new_node
         if (forward) then
            edge_length(attachment_edge) = attach_from_x - previous_sum
         else
            edge_length(attachment_edge) = path_sum - attach_from_x
         end if
         edge_count = edge_count + 1
         edge(edge_count, :) = [new_node, child]
         if (forward) then
            edge_length(edge_count) = path_sum - attach_from_x
         else
            edge_length(edge_count) = attach_from_x - previous_sum
         end if
         edge_count = edge_count + 1
         edge(edge_count, :) = [new_node, z]
         edge_length(edge_count) = minimum_distance
         k = k + 1
      end do

      tree = make_phylo_tree(n, edge(1:edge_count, :), edge_length(1:edge_count))
   end subroutine triang_mtds

   pure subroutine minimum_positive_pair(distance, x, y, info)
      !! Finds the first minimum positive off-diagonal distance using ape's scan order.
      real(dp), intent(in) :: distance(:, :) !! Complete square distance matrix.
      integer, intent(out) :: x !! First taxon index of the selected closest pair.
      integer, intent(out) :: y !! Second taxon index of the selected closest pair.
      integer, intent(out) :: info !! Zero on success; nonzero when no positive pair exists.
      integer :: i
      integer :: j
      integer :: n
      real(dp) :: minimum

      n = size(distance, 1)
      x = 0
      y = 0
      minimum = 0.0_dp
      do i = 1, n - 1
         do j = i + 1, n
            if (distance(i, j) <= 0.0_dp) cycle
            minimum = distance(i, j)
            x = i
            y = j
            exit
         end do
         if (x /= 0) exit
      end do
      if (x == 0) then
         info = 1
         return
      end if
      do i = 1, n - 1
         do j = i + 1, n
            if (distance(i, j) < minimum) then
               minimum = distance(i, j)
               x = i
               y = j
            end if
         end do
      end do
      info = 0
   end subroutine minimum_positive_pair

   pure subroutine partial_tree_path(x, y, max_node, edge, edge_count, path, info)
      !! Returns the node sequence between two vertices in a partially constructed rooted tree.
      integer, intent(in) :: x !! Starting vertex already present in the partial tree.
      integer, intent(in) :: y !! Ending vertex already present in the partial tree.
      integer, intent(in) :: max_node !! Largest node number addressable by the partial-tree storage.
      integer, intent(in) :: edge(:, :) !! Directed parent-child edges, with only the first `edge_count` rows active.
      integer, intent(in) :: edge_count !! Number of active rows in `edge`.
      integer, allocatable, intent(out) :: path(:) !! Node sequence from `x` through their common ancestor to `y`.
      integer, intent(out) :: info !! Zero on success; nonzero if either vertex is disconnected.
      integer, allocatable :: parent(:)
      integer, allocatable :: path_x(:)
      integer, allocatable :: path_y(:)
      integer :: i
      integer :: ix
      integer :: iy
      integer :: lca_x
      integer :: lca_y
      integer :: nx
      integer :: ny
      integer :: p

      allocate(parent(max_node), path_x(max_node), path_y(max_node))
      parent = 0
      do i = 1, edge_count
         if (edge(i, 2) >= 1 .and. edge(i, 2) <= max_node) parent(edge(i, 2)) = edge(i, 1)
      end do
      nx = 1
      path_x(nx) = x
      p = x
      do while (parent(p) /= 0)
         p = parent(p)
         nx = nx + 1
         path_x(nx) = p
      end do
      ny = 1
      path_y(ny) = y
      p = y
      do while (parent(p) /= 0)
         p = parent(p)
         ny = ny + 1
         path_y(ny) = p
      end do
      lca_x = 0
      lca_y = 0
      do ix = 1, nx
         do iy = 1, ny
            if (path_x(ix) == path_y(iy)) then
               lca_x = ix
               lca_y = iy
               exit
            end if
         end do
         if (lca_x /= 0) exit
      end do
      if (lca_x == 0) then
         allocate(path(0))
         info = 1
         return
      end if
      allocate(path(lca_x + lca_y - 1))
      path(1:lca_x) = path_x(1:lca_x)
      do i = 1, lca_y - 1
         path(lca_x + i) = path_y(lca_y - i)
      end do
      info = 0
   end subroutine partial_tree_path

   pure subroutine partial_edge_lookup(x, y, edge, edge_length, edge_count, index, length, forward)
      !! Locates an undirected partial-tree edge and reports traversal orientation.
      integer, intent(in) :: x !! First endpoint in the requested traversal direction.
      integer, intent(in) :: y !! Second endpoint in the requested traversal direction.
      integer, intent(in) :: edge(:, :) !! Directed parent-child edge storage.
      real(dp), intent(in) :: edge_length(:) !! Branch lengths aligned with `edge` rows.
      integer, intent(in) :: edge_count !! Number of active edge rows to search.
      integer, intent(out) :: index !! Matching edge row, or zero if the endpoints are not adjacent.
      real(dp), intent(out) :: length !! Length of the matching edge, or zero if not found.
      logical, intent(out) :: forward !! True when traversal is parent-to-child in the stored orientation.
      integer :: i

      index = 0
      length = 0.0_dp
      forward = .false.
      do i = 1, edge_count
         if (edge(i, 1) == x .and. edge(i, 2) == y) then
            index = i
            length = edge_length(i)
            forward = .true.
            return
         end if
         if (edge(i, 1) == y .and. edge(i, 2) == x) then
            index = i
            length = edge_length(i)
            forward = .false.
            return
         end if
      end do
   end subroutine partial_edge_lookup

   pure subroutine sort_by_known_count(known_count, order)
      !! Sorts taxa by decreasing known-distance counts with ape's stable pair-swap rule.
      integer, intent(inout) :: known_count(:) !! Per-taxon known-distance counts, reordered in place from largest to smallest.
      integer, intent(inout) :: order(:) !! Taxon labels permuted in lockstep with `known_count`.
      integer :: auxiliary
      integer :: i
      integer :: j

      do i = 1, size(known_count) - 1
         do j = i + 1, size(known_count)
            if (known_count(i) < known_count(j)) then
               auxiliary = known_count(i)
               known_count(i) = known_count(j)
               known_count(j) = auxiliary
               auxiliary = order(i)
               order(i) = order(j)
               order(j) = auxiliary
            end if
         end do
      end do
   end subroutine sort_by_known_count

   pure logical function valid_complete_distance_matrix(distance) result(ok)
      !! Validates a finite, symmetric, nonnegative complete distance matrix.
      real(dp), intent(in) :: distance(:, :) !! Candidate square taxon-distance matrix.
      integer :: i
      integer :: j
      integer :: n
      real(dp) :: tolerance

      n = size(distance, 1)
      ok = size(distance, 2) == n
      if (.not. ok) return
      tolerance = 64.0_dp * epsilon(1.0_dp)
      do i = 1, n - 1
         do j = i + 1, n
            if (ieee_is_nan(distance(i, j)) .or. ieee_is_nan(distance(j, i))) then
               ok = .false.
               return
            end if
            if (distance(i, j) < 0.0_dp .or. distance(j, i) < 0.0_dp) then
               ok = .false.
               return
            end if
            if (abs(distance(i, j) - distance(j, i)) > tolerance * &
               max(1.0_dp, abs(distance(i, j)), abs(distance(j, i)))) then
               ok = .false.
               return
            end if
         end do
      end do
   end function valid_complete_distance_matrix

   pure logical function valid_incomplete_distance_matrix_local(distance) result(ok)
      !! Validates symmetric known entries while permitting paired missing distances.
      real(dp), intent(in) :: distance(:, :) !! Square matrix; NaN or negative off-diagonal values mark missing pairs.
      integer :: i
      integer :: j
      integer :: n
      real(dp) :: tolerance

      n = size(distance, 1)
      ok = size(distance, 2) == n
      if (.not. ok) return
      tolerance = 64.0_dp * epsilon(1.0_dp)
      do i = 1, n - 1
         do j = i + 1, n
            if (missing_distance(distance(i, j)) .or. missing_distance(distance(j, i))) then
               if (.not. (missing_distance(distance(i, j)) .and. missing_distance(distance(j, i)))) then
                  ok = .false.
                  return
               end if
            else if (abs(distance(i, j) - distance(j, i)) > tolerance * &
               max(1.0_dp, abs(distance(i, j)), abs(distance(j, i)))) then
               ok = .false.
               return
            end if
         end do
      end do
   end function valid_incomplete_distance_matrix_local

   pure logical function all_distances_known(distance) result(known)
      !! Reports whether every off-diagonal distance is available.
      real(dp), intent(in) :: distance(:, :) !! Square matrix to scan for missing off-diagonal entries.
      integer :: i
      integer :: j

      known = .true.
      do i = 1, size(distance, 1) - 1
         do j = i + 1, size(distance, 1)
            if (missing_distance(distance(i, j))) then
               known = .false.
               return
            end if
         end do
      end do
   end function all_distances_known

   pure elemental logical function missing_distance(value) result(missing)
      !! Reports whether a triangles-method distance uses ape's missing-value convention.
      real(dp), intent(in) :: value !! Scalar distance; NaN or a negative value denotes missing data.

      missing = ieee_is_nan(value) .or. value < 0.0_dp
   end function missing_distance

end module ape_triangle_reconstruction
