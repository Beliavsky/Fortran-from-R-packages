! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
!
! FastME algorithms are derived from ape src/me.c, src/me_ols.c,
! src/me_balanced.c, src/NNI.c, src/bNNI.c, and src/SPR.c.
! Copyright 2007-2008 Olivier Gascuel, Rick Desper, Vincent Lefort,
! with modifications by Emmanuel Paradis; SPR code copyright 2009 Richard Desper.
!
! This modern Fortran implementation represents the evolving unrooted tree by
! an explicit edge graph and evaluates the same OLS/BME minimum-evolution
! objectives directly. It intentionally avoids the pointer-based C graph layer.
module ape_fastme
   use ieee_arithmetic, only : ieee_is_finite
   use r_kinds, only : dp
   use ape_types, only : phylo_tree, make_phylo_tree
   implicit none
   private

   real(dp), parameter :: fastme_epsilon = 1.0e-6_dp

   public :: fastme_ols
   public :: fastme_bal

contains

   pure subroutine fastme_ols(distance, tree, info, nni)
      !! Reconstructs a tree with ape's greedy OLS minimum-evolution algorithm and optional NNI search.
      real(dp), intent(in) :: distance(:, :) !! Symmetric taxon-distance matrix with at least three taxa.
      type(phylo_tree), intent(out) :: tree !! Reconstructed unrooted binary tree in ape-style directed encoding.
      integer, intent(out) :: info !! Status code: zero on success; nonzero for invalid distances or topology failure.
      logical, intent(in), optional :: nni !! Whether to run OLS nearest-neighbor interchange; defaults to true as in `fastme.ols`.
      integer, allocatable :: edges(:, :)
      logical :: do_nni
      integer :: n

      info = 0
      n = size(distance, 1)
      if (.not. valid_me_distance(distance)) then
         info = 1
         return
      end if

      do_nni = .true.
      if (present(nni)) do_nni = nni

      call greedy_insertion(distance, 1, edges, info)
      if (info /= 0) return
      if (do_nni) call nni_optimize(distance, 1, edges, info)
      if (info /= 0) return
      call graph_to_phylo(distance, edges, 1, tree, info)
   end subroutine fastme_ols

   pure subroutine fastme_bal(distance, tree, info, nni, spr)
      !! Reconstructs a tree with ape's balanced minimum-evolution insertion, bNNI, and optional SPR search.
      real(dp), intent(in) :: distance(:, :) !! Symmetric taxon-distance matrix with at least three taxa.
      type(phylo_tree), intent(out) :: tree !! Reconstructed unrooted binary tree in ape-style directed encoding.
      integer, intent(out) :: info !! Status code: zero on success; nonzero for invalid distances or topology failure.
      logical, intent(in), optional :: nni !! Whether to run balanced NNI optimization; defaults to true as in `fastme.bal`.
      logical, intent(in), optional :: spr !! Whether to run balanced SPR optimization; defaults to true as in `fastme.bal`.
      integer, allocatable :: edges(:, :)
      logical :: do_nni
      logical :: do_spr

      info = 0
      if (.not. valid_me_distance(distance)) then
         info = 1
         return
      end if

      do_nni = .true.
      if (present(nni)) do_nni = nni
      do_spr = .true.
      if (present(spr)) do_spr = spr

      call greedy_insertion(distance, 2, edges, info)
      if (info /= 0) return
      if (do_nni) call nni_optimize(distance, 2, edges, info)
      if (info /= 0) return
      if (do_spr) call spr_optimize(distance, edges, info)
      if (info /= 0) return
      call graph_to_phylo(distance, edges, 2, tree, info)
   end subroutine fastme_bal

   pure subroutine greedy_insertion(distance, method, edges, info)
      !! Builds the FastME starting topology by sequentially inserting taxa on the best current edge.
      real(dp), intent(in) :: distance(:, :) !! Complete symmetric distance matrix for all taxa.
      integer, intent(in) :: method !! Objective selector: 1 for OLS minimum evolution, 2 for balanced minimum evolution.
      integer, allocatable, intent(out) :: edges(:, :) !! Undirected binary-tree edge list after all taxa have been inserted.
      integer, intent(out) :: info !! Status code: zero on success; nonzero if objective evaluation fails.
      integer, allocatable :: candidate(:, :)
      integer, allocatable :: best_edges(:, :)
      integer :: edge_index
      integer :: k
      integer :: n
      real(dp) :: best_value
      real(dp) :: scale
      real(dp) :: tolerance
      real(dp) :: value

      info = 0
      n = size(distance, 1)
      allocate(edges(1, 2))
      edges(1, :) = [1, 2]
      call insert_tip(edges, 1, 3, n, candidate)
      call move_alloc(candidate, edges)

      do k = 4, n
         best_value = huge(1.0_dp)
         if (allocated(best_edges)) deallocate(best_edges)
         do edge_index = 1, size(edges, 1)
            call insert_tip(edges, edge_index, k, n, candidate)
            call me_objective(distance, candidate, k, method, value, info)
            if (info /= 0) return
            scale = max(1.0_dp, abs(best_value), abs(value))
            tolerance = 128.0_dp * epsilon(1.0_dp) * scale
            if (.not. allocated(best_edges) .or. value < best_value - tolerance) then
               best_value = value
               best_edges = candidate
            end if
         end do
         if (.not. allocated(best_edges)) then
            info = 2
            return
         end if
         call move_alloc(best_edges, edges)
      end do
   end subroutine greedy_insertion

   pure subroutine insert_tip(edges, edge_index, tip, n_total_tip, candidate)
      !! Splits one unrooted edge with a new internal node and attaches the next taxon as a pendant edge.
      integer, intent(in) :: edges(:, :) !! Current undirected edge list with shape `(2*k-3,2)` for the existing taxa.
      integer, intent(in) :: edge_index !! One-based row of `edges` to split.
      integer, intent(in) :: tip !! One-based taxon number being inserted; taxa are inserted in increasing order.
      integer, intent(in) :: n_total_tip !! Final number of taxa, used to reserve ape-compatible internal node numbers.
      integer, allocatable, intent(out) :: candidate(:, :) !! Edge list after insertion, containing two additional edges.
      integer :: i
      integer :: old_edges
      integer :: new_node
      integer :: out_index
      integer :: u
      integer :: v

      old_edges = size(edges, 1)
      new_node = n_total_tip + tip - 2
      u = edges(edge_index, 1)
      v = edges(edge_index, 2)
      allocate(candidate(old_edges + 2, 2))
      out_index = 0
      do i = 1, old_edges
         if (i /= edge_index) then
            out_index = out_index + 1
            candidate(out_index, :) = edges(i, :)
         else
            out_index = out_index + 1
            candidate(out_index, :) = [u, new_node]
            out_index = out_index + 1
            candidate(out_index, :) = [new_node, v]
            out_index = out_index + 1
            candidate(out_index, :) = [new_node, tip]
         end if
      end do
   end subroutine insert_tip

   pure subroutine nni_optimize(distance, method, edges, info)
      !! Repeatedly applies the best improving NNI move under the selected FastME objective.
      real(dp), intent(in) :: distance(:, :) !! Complete symmetric taxon-distance matrix.
      integer, intent(in) :: method !! Objective selector: 1 for OLS, 2 for balanced minimum evolution.
      integer, allocatable, intent(inout) :: edges(:, :) !! Undirected topology updated in place by accepted NNI moves.
      integer, intent(out) :: info !! Status code: zero on success; nonzero if an objective evaluation fails.
      integer, allocatable :: candidate(:, :)
      integer, allocatable :: best_edges(:, :)
      integer :: a_neighbors(3)
      integer :: b_neighbors(3)
      integer :: edge_index
      integer :: j
      integer :: n
      integer :: na
      integer :: nb
      integer :: swap_a
      integer :: swap_b
      integer :: u
      integer :: v
      real(dp) :: best_value
      real(dp) :: current_value
      real(dp) :: epsilon_search
      real(dp) :: scale
      real(dp) :: tolerance
      real(dp) :: value

      info = 0
      n = size(distance, 1)
      epsilon_search = sum(distance) / real(n * n, dp) * fastme_epsilon

      do
         call me_objective(distance, edges, n, method, current_value, info)
         if (info /= 0) return
         best_value = current_value
         if (allocated(best_edges)) deallocate(best_edges)

         do edge_index = 1, size(edges, 1)
            u = edges(edge_index, 1)
            v = edges(edge_index, 2)
            if (u <= n .or. v <= n) cycle
            call node_neighbors(edges, u, a_neighbors, na)
            call node_neighbors(edges, v, b_neighbors, nb)
            call remove_neighbor(a_neighbors, na, v)
            call remove_neighbor(b_neighbors, nb, u)
            if (na /= 2 .or. nb /= 2) cycle

            swap_a = a_neighbors(2)
            do j = 1, 2
               swap_b = b_neighbors(j)
               candidate = edges
               call swap_attachment(candidate, u, swap_a, v, swap_b, info)
               if (info /= 0) return
               call me_objective(distance, candidate, n, method, value, info)
               if (info /= 0) return
               scale = max(1.0_dp, abs(best_value), abs(value))
               tolerance = 128.0_dp * epsilon(1.0_dp) * scale
               if (value < best_value - tolerance) then
                  best_value = value
                  best_edges = candidate
               end if
            end do
         end do

         if (.not. allocated(best_edges)) exit
         if (best_value + epsilon_search >= current_value) exit
         call move_alloc(best_edges, edges)
      end do
   end subroutine nni_optimize

   pure subroutine spr_optimize(distance, edges, info)
      !! Repeatedly applies the best improving balanced-ME subtree-prune-and-regraft move.
      real(dp), intent(in) :: distance(:, :) !! Complete symmetric taxon-distance matrix.
      integer, allocatable, intent(inout) :: edges(:, :) !! Undirected topology updated in place by accepted SPR moves.
      integer, intent(out) :: info !! Status code: zero on success; nonzero if topology or objective evaluation fails.
      integer, allocatable :: base(:, :)
      integer, allocatable :: best_edges(:, :)
      integer, allocatable :: candidate(:, :)
      logical, allocatable :: pruned(:)
      integer :: cut_index
      integer :: direction
      integer :: n
      integer :: n_total
      integer :: a
      integer :: b
      integer :: neighbors(3)
      integer :: neighbor_count
      integer :: target_index
      integer :: x
      integer :: y
      real(dp) :: best_value
      real(dp) :: current_value
      real(dp) :: epsilon_search
      real(dp) :: scale
      real(dp) :: tolerance
      real(dp) :: value

      info = 0
      n = size(distance, 1)
      n_total = 2 * n - 2
      epsilon_search = sum(distance) / real(n * n, dp) * fastme_epsilon
      allocate(pruned(n_total))

      do
         call me_objective(distance, edges, n, 2, current_value, info)
         if (info /= 0) return
         best_value = current_value
         if (allocated(best_edges)) deallocate(best_edges)

         do cut_index = 1, size(edges, 1)
            do direction = 1, 2
               if (direction == 1) then
                  a = edges(cut_index, 1)
                  b = edges(cut_index, 2)
               else
                  a = edges(cut_index, 2)
                  b = edges(cut_index, 1)
               end if
               if (a <= n) cycle
               call node_neighbors(edges, a, neighbors, neighbor_count)
               call remove_neighbor(neighbors, neighbor_count, b)
               if (neighbor_count /= 2) cycle
               x = neighbors(1)
               y = neighbors(2)
               call component_nodes(edges, a, b, n_total, pruned)
               call suppress_cut_vertex(edges, a, b, x, y, base, info)
               if (info /= 0) return

               do target_index = 1, size(base, 1)
                  if (pruned(base(target_index, 1)) .or. pruned(base(target_index, 2))) cycle
                  call regraft_subtree(base, target_index, a, b, candidate)
                  call me_objective(distance, candidate, n, 2, value, info)
                  if (info /= 0) return
                  scale = max(1.0_dp, abs(best_value), abs(value))
                  tolerance = 128.0_dp * epsilon(1.0_dp) * scale
                  if (value < best_value - tolerance) then
                     best_value = value
                     best_edges = candidate
                  end if
               end do
            end do
         end do

         if (.not. allocated(best_edges)) exit
         if (best_value + epsilon_search >= current_value) exit
         call move_alloc(best_edges, edges)
      end do
   end subroutine spr_optimize

   pure subroutine suppress_cut_vertex(edges, a, b, x, y, base, info)
      !! Removes a prune edge and suppresses its degree-two retained-side endpoint before regrafting.
      integer, intent(in) :: edges(:, :) !! Current undirected binary-tree edge list.
      integer, intent(in) :: a !! Internal retained-side endpoint that becomes degree two after pruning.
      integer, intent(in) :: b !! Endpoint at the root of the component being pruned.
      integer, intent(in) :: x !! First retained neighbor of `a` other than `b`.
      integer, intent(in) :: y !! Second retained neighbor of `a` other than `b`.
      integer, allocatable, intent(out) :: base(:, :) !! Reduced retained/pruned forest plus the suppressing edge `x-y`.
      integer, intent(out) :: info !! Status code: zero on success; one if the three required edges were not found.
      integer :: i
      integer :: out_index
      integer :: removed

      info = 0
      allocate(base(size(edges, 1) - 2, 2))
      out_index = 0
      removed = 0
      do i = 1, size(edges, 1)
         if (same_edge(edges(i, 1), edges(i, 2), a, b) .or. &
            same_edge(edges(i, 1), edges(i, 2), a, x) .or. &
            same_edge(edges(i, 1), edges(i, 2), a, y)) then
            removed = removed + 1
            cycle
         end if
         out_index = out_index + 1
         base(out_index, :) = edges(i, :)
      end do
      if (removed /= 3) then
         info = 1
         return
      end if
      out_index = out_index + 1
      base(out_index, :) = [x, y]
   end subroutine suppress_cut_vertex

   pure subroutine regraft_subtree(base, target_index, attach_node, pruned_root, candidate)
      !! Subdivides one retained-tree edge with the suppressed node and reconnects the pruned component.
      integer, intent(in) :: base(:, :) !! Forest edge list after pruning and suppression.
      integer, intent(in) :: target_index !! Row of `base` whose retained edge is subdivided for regrafting.
      integer, intent(in) :: attach_node !! Reused internal node that subdivides the target edge.
      integer, intent(in) :: pruned_root !! Node at the root of the pruned component to reconnect.
      integer, allocatable, intent(out) :: candidate(:, :) !! Full binary-tree edge list after the proposed regraft.
      integer :: c
      integer :: d
      integer :: i
      integer :: out_index

      c = base(target_index, 1)
      d = base(target_index, 2)
      allocate(candidate(size(base, 1) + 2, 2))
      out_index = 0
      do i = 1, size(base, 1)
         if (i /= target_index) then
            out_index = out_index + 1
            candidate(out_index, :) = base(i, :)
         else
            out_index = out_index + 1
            candidate(out_index, :) = [c, attach_node]
            out_index = out_index + 1
            candidate(out_index, :) = [attach_node, d]
            out_index = out_index + 1
            candidate(out_index, :) = [attach_node, pruned_root]
         end if
      end do
   end subroutine regraft_subtree

   pure subroutine swap_attachment(edges, u, a, v, b, info)
      !! Performs one NNI alternative by exchanging the attachments `u-a` and `v-b` across internal edge `u-v`.
      integer, intent(inout) :: edges(:, :) !! Undirected topology modified in place.
      integer, intent(in) :: u !! First endpoint of the central internal edge.
      integer, intent(in) :: a !! Neighbor detached from `u` and attached to `v`.
      integer, intent(in) :: v !! Second endpoint of the central internal edge.
      integer, intent(in) :: b !! Neighbor detached from `v` and attached to `u`.
      integer, intent(out) :: info !! Status code: zero on success; one if either attachment edge is absent.
      integer :: edge_a
      integer :: edge_b

      info = 0
      edge_a = find_edge(edges, u, a)
      edge_b = find_edge(edges, v, b)
      if (edge_a == 0 .or. edge_b == 0) then
         info = 1
         return
      end if
      call replace_edge_endpoint(edges, edge_a, u, v)
      call replace_edge_endpoint(edges, edge_b, v, u)
   end subroutine swap_attachment

   pure subroutine replace_edge_endpoint(edges, edge_index, old_node, new_node)
      !! Replaces one endpoint of an undirected edge while preserving row orientation when possible.
      integer, intent(inout) :: edges(:, :) !! Undirected edge matrix containing the row to edit.
      integer, intent(in) :: edge_index !! One-based row of `edges` whose endpoint is replaced.
      integer, intent(in) :: old_node !! Endpoint value that must currently occur in the selected row.
      integer, intent(in) :: new_node !! Replacement node number for that endpoint.

      if (edges(edge_index, 1) == old_node) then
         edges(edge_index, 1) = new_node
      else if (edges(edge_index, 2) == old_node) then
         edges(edge_index, 2) = new_node
      end if
   end subroutine replace_edge_endpoint

   pure subroutine me_objective(distance, edges, n_tip, method, value, info)
      !! Returns total fitted branch length for an OLS or balanced-ME topology.
      real(dp), intent(in) :: distance(:, :) !! Full distance matrix whose leading `n_tip` taxa are active.
      integer, intent(in) :: edges(:, :) !! Undirected edge list for the active topology.
      integer, intent(in) :: n_tip !! Number of taxa currently represented by `edges`.
      integer, intent(in) :: method !! Objective selector: 1 for OLS, 2 for balanced minimum evolution.
      real(dp), intent(out) :: value !! Sum of fitted branch lengths for the requested topology.
      integer, intent(out) :: info !! Status code: zero on success; nonzero for malformed topology or unknown method.
      real(dp), allocatable :: lengths(:)

      value = huge(1.0_dp)
      select case (method)
      case (1)
         call ols_branch_lengths(distance, edges, n_tip, lengths, info)
      case (2)
         call balanced_branch_lengths(distance, edges, n_tip, lengths, info)
      case default
         info = 1
         return
      end select
      if (info /= 0) return
      value = sum(lengths)
   end subroutine me_objective

   pure subroutine ols_branch_lengths(distance, edges, n_tip, lengths, info)
      !! Computes FastME OLS edge lengths from arithmetic subtree-average distances.
      real(dp), intent(in) :: distance(:, :) !! Full taxon-distance matrix whose leading `n_tip` rows/columns are active.
      integer, intent(in) :: edges(:, :) !! Undirected binary-tree edge list.
      integer, intent(in) :: n_tip !! Number of active taxa represented by the topology.
      real(dp), allocatable, intent(out) :: lengths(:) !! OLS length for each row of `edges`, in matching order.
      integer, intent(out) :: info !! Status code: zero on success; nonzero if the edge graph is not binary.
      logical, allocatable :: mask_a(:)
      logical, allocatable :: mask_b(:)
      logical, allocatable :: mask_c(:)
      logical, allocatable :: mask_d(:)
      integer :: count_a
      integer :: count_b
      integer :: count_c
      integer :: count_d
      integer :: i
      integer :: neighbors_u(3)
      integer :: neighbors_v(3)
      integer :: nu
      integer :: nv
      integer :: tip
      integer :: u
      integer :: v
      integer :: x
      real(dp) :: d_ab
      real(dp) :: d_ac
      real(dp) :: d_ad
      real(dp) :: d_bc
      real(dp) :: d_bd
      real(dp) :: d_cd
      real(dp) :: lambda

      info = 0
      allocate(lengths(size(edges, 1)))
      allocate(mask_a(n_tip), mask_b(n_tip), mask_c(n_tip), mask_d(n_tip))
      do i = 1, size(edges, 1)
         u = edges(i, 1)
         v = edges(i, 2)
         if (u <= n_tip .or. v <= n_tip) then
            if (u <= n_tip) then
               tip = u
               x = v
            else
               tip = v
               x = u
            end if
            call node_neighbors(edges, x, neighbors_u, nu)
            call remove_neighbor(neighbors_u, nu, tip)
            if (nu /= 2) then
               info = 1
               return
            end if
            call component_tip_mask(edges, x, neighbors_u(1), n_tip, mask_a)
            call component_tip_mask(edges, x, neighbors_u(2), n_tip, mask_b)
            mask_c = .false.
            mask_c(tip) = .true.
            d_ac = uniform_average(distance, mask_c, mask_a)
            d_bc = uniform_average(distance, mask_c, mask_b)
            d_ab = uniform_average(distance, mask_a, mask_b)
            lengths(i) = 0.5_dp * (d_ac + d_bc - d_ab)
         else
            call node_neighbors(edges, u, neighbors_u, nu)
            call node_neighbors(edges, v, neighbors_v, nv)
            call remove_neighbor(neighbors_u, nu, v)
            call remove_neighbor(neighbors_v, nv, u)
            if (nu /= 2 .or. nv /= 2) then
               info = 1
               return
            end if
            call component_tip_mask(edges, u, neighbors_u(1), n_tip, mask_a)
            call component_tip_mask(edges, u, neighbors_u(2), n_tip, mask_b)
            call component_tip_mask(edges, v, neighbors_v(1), n_tip, mask_c)
            call component_tip_mask(edges, v, neighbors_v(2), n_tip, mask_d)
            count_a = count(mask_a)
            count_b = count(mask_b)
            count_c = count(mask_c)
            count_d = count(mask_d)
            if (count_a == 0 .or. count_b == 0 .or. count_c == 0 .or. count_d == 0) then
               info = 1
               return
            end if
            lambda = real(count_a * count_c + count_d * count_b, dp) / &
               real((count_c + count_d) * (count_a + count_b), dp)
            d_ab = uniform_average(distance, mask_a, mask_b)
            d_ac = uniform_average(distance, mask_a, mask_c)
            d_ad = uniform_average(distance, mask_a, mask_d)
            d_bc = uniform_average(distance, mask_b, mask_c)
            d_bd = uniform_average(distance, mask_b, mask_d)
            d_cd = uniform_average(distance, mask_c, mask_d)
            lengths(i) = 0.5_dp * (lambda * (d_bc + d_ad) + (1.0_dp - lambda) * (d_ac + d_bd) - &
               (d_cd + d_ab))
         end if
      end do
   end subroutine ols_branch_lengths

   pure subroutine balanced_branch_lengths(distance, edges, n_tip, lengths, info)
      !! Computes FastME balanced edge lengths from recursively half-weighted subtree averages.
      real(dp), intent(in) :: distance(:, :) !! Full taxon-distance matrix whose leading `n_tip` rows/columns are active.
      integer, intent(in) :: edges(:, :) !! Undirected binary-tree edge list.
      integer, intent(in) :: n_tip !! Number of active taxa represented by the topology.
      real(dp), allocatable, intent(out) :: lengths(:) !! Balanced-ME length for each row of `edges`, in matching order.
      integer, intent(out) :: info !! Status code: zero on success; nonzero if the edge graph is not binary.
      real(dp), allocatable :: weight_a(:)
      real(dp), allocatable :: weight_b(:)
      real(dp), allocatable :: weight_c(:)
      real(dp), allocatable :: weight_d(:)
      integer :: i
      integer :: neighbors_u(3)
      integer :: neighbors_v(3)
      integer :: nu
      integer :: nv
      integer :: tip
      integer :: u
      integer :: v
      integer :: x
      real(dp) :: d_ab
      real(dp) :: d_ac
      real(dp) :: d_ad
      real(dp) :: d_bc
      real(dp) :: d_bd
      real(dp) :: d_cd

      info = 0
      allocate(lengths(size(edges, 1)))
      allocate(weight_a(n_tip), weight_b(n_tip), weight_c(n_tip), weight_d(n_tip))
      do i = 1, size(edges, 1)
         u = edges(i, 1)
         v = edges(i, 2)
         if (u <= n_tip .or. v <= n_tip) then
            if (u <= n_tip) then
               tip = u
               x = v
            else
               tip = v
               x = u
            end if
            call node_neighbors(edges, x, neighbors_u, nu)
            call remove_neighbor(neighbors_u, nu, tip)
            if (nu /= 2) then
               info = 1
               return
            end if
            call component_balanced_weights(edges, x, neighbors_u(1), n_tip, weight_a, info)
            if (info /= 0) return
            call component_balanced_weights(edges, x, neighbors_u(2), n_tip, weight_b, info)
            if (info /= 0) return
            weight_c = 0.0_dp
            weight_c(tip) = 1.0_dp
            d_ac = weighted_average(distance, weight_c, weight_a)
            d_bc = weighted_average(distance, weight_c, weight_b)
            d_ab = weighted_average(distance, weight_a, weight_b)
            lengths(i) = 0.5_dp * (d_ac + d_bc - d_ab)
         else
            call node_neighbors(edges, u, neighbors_u, nu)
            call node_neighbors(edges, v, neighbors_v, nv)
            call remove_neighbor(neighbors_u, nu, v)
            call remove_neighbor(neighbors_v, nv, u)
            if (nu /= 2 .or. nv /= 2) then
               info = 1
               return
            end if
            call component_balanced_weights(edges, u, neighbors_u(1), n_tip, weight_a, info)
            if (info /= 0) return
            call component_balanced_weights(edges, u, neighbors_u(2), n_tip, weight_b, info)
            if (info /= 0) return
            call component_balanced_weights(edges, v, neighbors_v(1), n_tip, weight_c, info)
            if (info /= 0) return
            call component_balanced_weights(edges, v, neighbors_v(2), n_tip, weight_d, info)
            if (info /= 0) return
            d_ab = weighted_average(distance, weight_a, weight_b)
            d_ac = weighted_average(distance, weight_a, weight_c)
            d_ad = weighted_average(distance, weight_a, weight_d)
            d_bc = weighted_average(distance, weight_b, weight_c)
            d_bd = weighted_average(distance, weight_b, weight_d)
            d_cd = weighted_average(distance, weight_c, weight_d)
            lengths(i) = 0.25_dp * (d_ac + d_ad + d_bc + d_bd) - 0.5_dp * (d_ab + d_cd)
         end if
      end do
   end subroutine balanced_branch_lengths

   pure subroutine component_tip_mask(edges, previous, start, n_tip, mask)
      !! Marks the tips in the connected component entered through one oriented side of an undirected edge.
      integer, intent(in) :: edges(:, :) !! Undirected tree edge list.
      integer, intent(in) :: previous !! Boundary node excluded from traversal so the selected edge is effectively cut.
      integer, intent(in) :: start !! First node on the component side to traverse.
      integer, intent(in) :: n_tip !! Number of terminal taxa; tip node numbers are `1:n_tip`.
      logical, intent(out) :: mask(:) !! Logical tip-membership vector of length `n_tip`.
      integer, allocatable :: parent(:)
      integer, allocatable :: queue(:)
      integer :: current
      integer :: head
      integer :: i
      integer :: n_total
      integer :: tail
      integer :: next

      n_total = max(maxval(edges), n_tip)
      allocate(parent(n_total), queue(n_total))
      parent = 0
      mask = .false.
      head = 1
      tail = 1
      queue(1) = start
      parent(start) = previous
      do while (head <= tail)
         current = queue(head)
         head = head + 1
         if (current <= n_tip) mask(current) = .true.
         do i = 1, size(edges, 1)
            next = 0
            if (edges(i, 1) == current) then
               next = edges(i, 2)
            else if (edges(i, 2) == current) then
               next = edges(i, 1)
            end if
            if (next == 0 .or. next == parent(current)) cycle
            if (parent(next) /= 0 .or. next == start) cycle
            parent(next) = current
            tail = tail + 1
            queue(tail) = next
         end do
      end do
   end subroutine component_tip_mask

   pure subroutine component_balanced_weights(edges, previous, start, n_tip, weights, info)
      !! Computes FastME's recursively half-weighted distribution over tips on one oriented component.
      integer, intent(in) :: edges(:, :) !! Undirected binary-tree edge list.
      integer, intent(in) :: previous !! Boundary node excluded from traversal so the selected edge is effectively cut.
      integer, intent(in) :: start !! First node on the component side to traverse.
      integer, intent(in) :: n_tip !! Number of terminal taxa; tip node numbers are `1:n_tip`.
      real(dp), intent(out) :: weights(:) !! Tip weights of length `n_tip`, summing to one for a valid binary component.
      integer, intent(out) :: info !! Status code: zero on success; nonzero if an internal component node does not bifurcate.
      integer, allocatable :: parent_stack(:)
      integer, allocatable :: node_stack(:)
      real(dp), allocatable :: coefficient_stack(:)
      integer :: current
      integer :: i
      integer :: neighbors(3)
      integer :: neighbor_count
      integer :: parent
      integer :: stack_size
      real(dp) :: coefficient

      info = 0
      weights = 0.0_dp
      allocate(parent_stack(max(1, 2 * n_tip - 2)))
      allocate(node_stack(max(1, 2 * n_tip - 2)))
      allocate(coefficient_stack(max(1, 2 * n_tip - 2)))
      stack_size = 1
      parent_stack(1) = previous
      node_stack(1) = start
      coefficient_stack(1) = 1.0_dp

      do while (stack_size > 0)
         parent = parent_stack(stack_size)
         current = node_stack(stack_size)
         coefficient = coefficient_stack(stack_size)
         stack_size = stack_size - 1
         if (current <= n_tip) then
            weights(current) = weights(current) + coefficient
            cycle
         end if
         call node_neighbors(edges, current, neighbors, neighbor_count)
         call remove_neighbor(neighbors, neighbor_count, parent)
         if (neighbor_count /= 2) then
            info = 1
            return
         end if
         do i = 1, 2
            stack_size = stack_size + 1
            parent_stack(stack_size) = current
            node_stack(stack_size) = neighbors(i)
            coefficient_stack(stack_size) = 0.5_dp * coefficient
         end do
      end do
   end subroutine component_balanced_weights

   pure real(dp) function uniform_average(distance, mask_a, mask_b) result(value)
      !! Returns the arithmetic mean of all pairwise distances between two nonempty tip sets.
      real(dp), intent(in) :: distance(:, :) !! Taxon-distance matrix containing both tip sets.
      logical, intent(in) :: mask_a(:) !! Membership mask for the first tip set.
      logical, intent(in) :: mask_b(:) !! Membership mask for the second tip set.
      integer :: i
      integer :: j
      integer :: pairs

      value = 0.0_dp
      pairs = 0
      do i = 1, size(mask_a)
         if (.not. mask_a(i)) cycle
         do j = 1, size(mask_b)
            if (.not. mask_b(j)) cycle
            value = value + distance(i, j)
            pairs = pairs + 1
         end do
      end do
      if (pairs > 0) value = value / real(pairs, dp)
   end function uniform_average

   pure real(dp) function weighted_average(distance, weight_a, weight_b) result(value)
      !! Returns a pairwise distance average under two normalized FastME subtree-weight distributions.
      real(dp), intent(in) :: distance(:, :) !! Taxon-distance matrix containing both weighted tip sets.
      real(dp), intent(in) :: weight_a(:) !! Nonnegative weights for the first component, normally summing to one.
      real(dp), intent(in) :: weight_b(:) !! Nonnegative weights for the second component, normally summing to one.
      integer :: i
      integer :: j

      value = 0.0_dp
      do i = 1, size(weight_a)
         if (abs(weight_a(i)) <= tiny(1.0_dp)) cycle
         do j = 1, size(weight_b)
            if (abs(weight_b(j)) <= tiny(1.0_dp)) cycle
            value = value + weight_a(i) * weight_b(j) * distance(i, j)
         end do
      end do
   end function weighted_average

   pure subroutine graph_to_phylo(distance, edges, method, tree, info)
      !! Orients a final unrooted FastME graph from internal node `n_tip+1` and attaches fitted branch lengths.
      real(dp), intent(in) :: distance(:, :) !! Complete taxon-distance matrix used to fit final branch lengths.
      integer, intent(in) :: edges(:, :) !! Final undirected binary-tree edge list.
      integer, intent(in) :: method !! Length selector: 1 for OLS or 2 for balanced minimum evolution.
      type(phylo_tree), intent(out) :: tree !! Ape-style directed tree rooted at the first internal node.
      integer, intent(out) :: info !! Status code: zero on success; nonzero for malformed topology or unknown method.
      integer, allocatable :: directed(:, :)
      integer, allocatable :: parent(:)
      integer, allocatable :: queue(:)
      real(dp), allocatable :: lengths(:)
      integer :: current
      integer :: head
      integer :: i
      integer :: n
      integer :: n_total
      integer :: next
      integer :: root
      integer :: tail

      info = 0
      n = size(distance, 1)
      n_total = 2 * n - 2
      root = n + 1
      select case (method)
      case (1)
         call ols_branch_lengths(distance, edges, n, lengths, info)
      case (2)
         call balanced_branch_lengths(distance, edges, n, lengths, info)
      case default
         info = 1
         return
      end select
      if (info /= 0) return

      allocate(parent(n_total), queue(n_total), directed(size(edges, 1), 2))
      parent = 0
      head = 1
      tail = 1
      queue(1) = root
      parent(root) = -1
      do while (head <= tail)
         current = queue(head)
         head = head + 1
         do i = 1, size(edges, 1)
            next = 0
            if (edges(i, 1) == current) then
               next = edges(i, 2)
            else if (edges(i, 2) == current) then
               next = edges(i, 1)
            end if
            if (next == 0) cycle
            if (parent(next) /= 0) cycle
            parent(next) = current
            tail = tail + 1
            queue(tail) = next
         end do
      end do
      if (tail /= n_total) then
         info = 2
         return
      end if

      do i = 1, size(edges, 1)
         if (parent(edges(i, 2)) == edges(i, 1)) then
            directed(i, :) = edges(i, :)
         else if (parent(edges(i, 1)) == edges(i, 2)) then
            directed(i, :) = [edges(i, 2), edges(i, 1)]
         else
            info = 2
            return
         end if
      end do
      tree = make_phylo_tree(n, directed, lengths)
   end subroutine graph_to_phylo

   pure subroutine component_nodes(edges, previous, start, n_total, membership)
      !! Marks every node in one component obtained by cutting the oriented edge `previous-start`.
      integer, intent(in) :: edges(:, :) !! Undirected tree edge list.
      integer, intent(in) :: previous !! Boundary node on the excluded side of the cut edge.
      integer, intent(in) :: start !! First node in the component to mark.
      integer, intent(in) :: n_total !! Total number of tip and internal nodes in the final binary tree.
      logical, intent(out) :: membership(:) !! Node-membership vector of length `n_total` for the selected component.
      integer, allocatable :: parent(:)
      integer, allocatable :: queue(:)
      integer :: current
      integer :: head
      integer :: i
      integer :: next
      integer :: tail

      membership = .false.
      allocate(parent(n_total), queue(n_total))
      parent = 0
      head = 1
      tail = 1
      queue(1) = start
      parent(start) = previous
      membership(start) = .true.
      do while (head <= tail)
         current = queue(head)
         head = head + 1
         do i = 1, size(edges, 1)
            next = 0
            if (edges(i, 1) == current) then
               next = edges(i, 2)
            else if (edges(i, 2) == current) then
               next = edges(i, 1)
            end if
            if (next == 0 .or. next == parent(current)) cycle
            if (membership(next)) cycle
            parent(next) = current
            membership(next) = .true.
            tail = tail + 1
            queue(tail) = next
         end do
      end do
   end subroutine component_nodes

   pure subroutine node_neighbors(edges, node, neighbors, count_neighbors)
      !! Lists up to three neighbors of a node in edge-row order for an unrooted binary topology.
      integer, intent(in) :: edges(:, :) !! Undirected binary-tree edge list.
      integer, intent(in) :: node !! Node number whose incident neighbors are requested.
      integer, intent(out) :: neighbors(3) !! Neighbor node numbers in the order encountered in `edges`.
      integer, intent(out) :: count_neighbors !! Number of incident neighbors found, between zero and three for valid input.
      integer :: i

      neighbors = 0
      count_neighbors = 0
      do i = 1, size(edges, 1)
         if (edges(i, 1) == node) then
            count_neighbors = count_neighbors + 1
            if (count_neighbors <= 3) neighbors(count_neighbors) = edges(i, 2)
         else if (edges(i, 2) == node) then
            count_neighbors = count_neighbors + 1
            if (count_neighbors <= 3) neighbors(count_neighbors) = edges(i, 1)
         end if
      end do
   end subroutine node_neighbors

   pure subroutine remove_neighbor(neighbors, count_neighbors, excluded)
      !! Removes one specified node from a short neighbor list while preserving the relative order of the others.
      integer, intent(inout) :: neighbors(3) !! Neighbor list to edit in place.
      integer, intent(inout) :: count_neighbors !! Number of valid entries in `neighbors`, decremented if `excluded` is found.
      integer, intent(in) :: excluded !! Node number to remove from the valid prefix of `neighbors`.
      integer :: i
      integer :: j

      do i = 1, count_neighbors
         if (neighbors(i) /= excluded) cycle
         do j = i, count_neighbors - 1
            neighbors(j) = neighbors(j + 1)
         end do
         neighbors(count_neighbors) = 0
         count_neighbors = count_neighbors - 1
         return
      end do
   end subroutine remove_neighbor

   pure integer function find_edge(edges, u, v) result(index)
      !! Finds the row containing an undirected node pair, returning zero if no such edge exists.
      integer, intent(in) :: edges(:, :) !! Undirected edge list to search.
      integer, intent(in) :: u !! First endpoint of the requested edge.
      integer, intent(in) :: v !! Second endpoint of the requested edge.
      integer :: i

      index = 0
      do i = 1, size(edges, 1)
         if (same_edge(edges(i, 1), edges(i, 2), u, v)) then
            index = i
            return
         end if
      end do
   end function find_edge

   pure logical function same_edge(edge_u, edge_v, u, v) result(matches)
      !! Tests two stored endpoints against an undirected node pair.
      integer, intent(in) :: edge_u !! First stored endpoint of the edge.
      integer, intent(in) :: edge_v !! Second stored endpoint of the edge.
      integer, intent(in) :: u !! First endpoint of the requested undirected pair.
      integer, intent(in) :: v !! Second endpoint of the requested undirected pair.

      matches = (edge_u == u .and. edge_v == v) .or. (edge_u == v .and. edge_v == u)
   end function same_edge

   pure logical function valid_me_distance(distance) result(ok)
      !! Checks FastME's matrix shape, finite values, zero diagonal, nonnegative distances, and symmetry.
      real(dp), intent(in) :: distance(:, :) !! Candidate square taxon-distance matrix.
      integer :: i
      integer :: j
      integer :: n
      real(dp) :: scale

      ok = .false.
      n = size(distance, 1)
      if (n < 3 .or. size(distance, 2) /= n) return
      do i = 1, n
         if (.not. ieee_is_finite(distance(i, i))) return
         if (abs(distance(i, i)) > 128.0_dp * epsilon(1.0_dp)) return
         do j = i + 1, n
            if (.not. ieee_is_finite(distance(i, j)) .or. .not. ieee_is_finite(distance(j, i))) return
            if (distance(i, j) < 0.0_dp .or. distance(j, i) < 0.0_dp) return
            scale = max(1.0_dp, abs(distance(i, j)), abs(distance(j, i)))
            if (abs(distance(i, j) - distance(j, i)) > 128.0_dp * epsilon(1.0_dp) * scale) return
         end do
      end do
      ok = .true.
   end function valid_me_distance

end module ape_fastme
