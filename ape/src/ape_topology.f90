! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
!
! Topology utilities translate computational ideas from ape R/is.binary.tree.R,
! R/is.ultrametric.R, R/vcv.phylo.R, R/is.monophyletic.R, R/which.edge.R,
! R/dist.topo.R (Copyright 2005-2023 Emmanuel Paradis, 2016-2021 Klaus
! Schliep), and src/prop_part.cpp (Copyright 2017 Klaus Schliep).
module ape_topology
   use r_kinds, only : dp
   use ape_types, only : phylo_tree, child_counts, parent_vector, edge_index_to_child
   use ape_tree_algorithms, only : node_depth_edgelength, mrca
   implicit none
   private

   public :: is_binary_tree
   public :: is_ultrametric_tree
   public :: phylogenetic_vcv
   public :: mrca_many
   public :: is_monophyletic
   public :: tip_descendant_matrix
   public :: clade_tips
   public :: topological_distance_ph85
   public :: branch_score_distance

contains

   pure elemental logical function is_binary_tree(tree, allow_trifurcating_root) result(binary)
      !! Tests whether every internal node is dichotomous, optionally allowing an unrooted trifurcating root.
      type(phylo_tree), intent(in) :: tree !! Tree whose internal-node out-degrees are inspected.
      logical, intent(in), optional :: allow_trifurcating_root !! Accept a three-child root as an unrooted binary tree.
      integer, allocatable :: counts(:)
      integer :: node
      integer :: root_node
      logical :: allow_three

      binary = .false.
      if (.not. tree%valid()) return
      counts = child_counts(tree)
      root_node = tree%root()
      allow_three = .false.
      if (present(allow_trifurcating_root)) allow_three = allow_trifurcating_root
      do node = tree%n_tip + 1, tree%total_nodes()
         if (node == root_node .and. allow_three) then
            if (counts(node) /= 2 .and. counts(node) /= 3) return
         else
            if (counts(node) /= 2) return
         end if
      end do
      binary = .true.
   end function is_binary_tree

   pure elemental logical function is_ultrametric_tree(tree, tolerance, option) result(ultrametric)
      !! Tests equality of root-to-tip distances using ape's relative-range or variance criterion.
      type(phylo_tree), intent(in) :: tree !! Rooted tree with branch lengths.
      real(dp), intent(in), optional :: tolerance !! Nonnegative tolerance; default is `sqrt(epsilon(1.0_dp))`.
      integer, intent(in), optional :: option !! Criterion: 1 is relative range, 2 is sample variance.
      real(dp), allocatable :: depth(:)
      real(dp) :: criterion
      real(dp) :: mean_depth
      real(dp) :: tol
      real(dp) :: maximum
      real(dp) :: minimum
      integer :: info
      integer :: method

      ultrametric = .false.
      call node_depth_edgelength(tree, depth, info)
      if (info /= 0 .or. tree%n_tip < 1) return
      tol = sqrt(epsilon(1.0_dp))
      if (present(tolerance)) tol = tolerance
      if (tol < 0.0_dp) return
      method = 1
      if (present(option)) method = option
      select case (method)
      case (1)
         minimum = minval(depth(1:tree%n_tip))
         maximum = maxval(depth(1:tree%n_tip))
         if (abs(maximum) <= tiny(1.0_dp)) then
            criterion = maximum - minimum
         else
            criterion = (maximum - minimum) / abs(maximum)
         end if
      case (2)
         if (tree%n_tip == 1) then
            criterion = 0.0_dp
         else
            mean_depth = sum(depth(1:tree%n_tip)) / real(tree%n_tip, dp)
            criterion = sum((depth(1:tree%n_tip) - mean_depth)**2) / real(tree%n_tip - 1, dp)
         end if
      case default
         return
      end select
      ultrametric = abs(criterion) <= tol
   end function is_ultrametric_tree

   pure subroutine phylogenetic_vcv(tree, covariance, info, correlation)
      !! Computes the Brownian-motion phylogenetic variance-covariance or correlation matrix for the tips.
      type(phylo_tree), intent(in) :: tree !! Rooted tree with branch lengths.
      real(dp), allocatable, intent(out) :: covariance(:, :) !! Tip covariance matrix whose covariance is root-to-MRCA path length.
      integer, intent(out) :: info !! Status code: zero on success, nonzero if tree depths or correlations are undefined.
      logical, intent(in), optional :: correlation !! If true, normalize the covariance matrix to unit diagonal.
      real(dp), allocatable :: depth(:)
      real(dp), allocatable :: inverse_sd(:)
      integer :: ancestor
      integer :: i
      integer :: j
      logical :: as_correlation

      call node_depth_edgelength(tree, depth, info)
      allocate(covariance(tree%n_tip, tree%n_tip))
      covariance = 0.0_dp
      if (info /= 0) return
      do i = 1, tree%n_tip
         covariance(i, i) = depth(i)
         do j = i + 1, tree%n_tip
            ancestor = mrca(tree, i, j)
            if (ancestor == 0) then
               info = 2
               return
            end if
            covariance(i, j) = depth(ancestor)
            covariance(j, i) = covariance(i, j)
         end do
      end do
      as_correlation = .false.
      if (present(correlation)) as_correlation = correlation
      if (.not. as_correlation) return
      allocate(inverse_sd(tree%n_tip))
      do i = 1, tree%n_tip
         if (covariance(i, i) <= 0.0_dp) then
            info = 3
            return
         end if
         inverse_sd(i) = 1.0_dp / sqrt(covariance(i, i))
      end do
      do i = 1, tree%n_tip
         do j = 1, tree%n_tip
            covariance(i, j) = inverse_sd(i) * covariance(i, j) * inverse_sd(j)
         end do
         covariance(i, i) = 1.0_dp
      end do
   end subroutine phylogenetic_vcv

   pure integer function mrca_many(tree, nodes) result(ancestor)
      !! Finds the most recent common ancestor shared by every node in a nonempty set.
      type(phylo_tree), intent(in) :: tree !! Rooted tree containing the requested nodes.
      integer, intent(in) :: nodes(:) !! Ape node numbers whose joint most recent common ancestor is requested.
      integer :: i

      ancestor = 0
      if (size(nodes) == 0) return
      if (.not. tree%valid()) return
      ancestor = nodes(1)
      if (ancestor < 1 .or. ancestor > tree%total_nodes()) then
         ancestor = 0
         return
      end if
      do i = 2, size(nodes)
         ancestor = mrca(tree, ancestor, nodes(i))
         if (ancestor == 0) return
      end do
   end function mrca_many

   pure logical function is_monophyletic(tree, tips) result(monophyletic)
      !! Tests whether a set of tip numbers is exactly the descendant-tip set of its MRCA.
      type(phylo_tree), intent(in) :: tree !! Rooted tree in ape node numbering.
      integer, intent(in) :: tips(:) !! Unique tip numbers in the range `1:n_tip`.
      logical, allocatable :: wanted(:)
      logical, allocatable :: descendants(:)
      integer :: ancestor
      integer :: i

      monophyletic = .false.
      if (.not. tree%valid()) return
      if (size(tips) == 0) return
      allocate(wanted(tree%n_tip))
      wanted = .false.
      do i = 1, size(tips)
         if (tips(i) < 1 .or. tips(i) > tree%n_tip) return
         if (wanted(tips(i))) return
         wanted(tips(i)) = .true.
      end do
      if (size(tips) == 1 .or. size(tips) == tree%n_tip) then
         monophyletic = .true.
         return
      end if
      ancestor = mrca_many(tree, tips)
      if (ancestor == 0) return
      call clade_tips(tree, ancestor, descendants)
      monophyletic = all(descendants .eqv. wanted)
   end function is_monophyletic

   pure subroutine tip_descendant_matrix(tree, descendants, info)
      !! Builds a logical matrix identifying the tip descendants of every internal node.
      type(phylo_tree), intent(in) :: tree !! Rooted tree whose internal clades are encoded.
      logical, allocatable, intent(out) :: descendants(:, :) !! Shape `(n_node, n_tip)` with one internal clade per row.
      integer, intent(out) :: info !! Status code: zero on success, nonzero for an invalid tree.
      integer, allocatable :: parent(:)
      integer :: current
      integer :: tip

      allocate(descendants(tree%n_node, tree%n_tip))
      descendants = .false.
      info = 0
      if (.not. tree%valid()) then
         info = 1
         return
      end if
      parent = parent_vector(tree)
      do tip = 1, tree%n_tip
         current = parent(tip)
         do while (current /= 0)
            descendants(current - tree%n_tip, tip) = .true.
            current = parent(current)
         end do
      end do
   end subroutine tip_descendant_matrix

   pure subroutine clade_tips(tree, node, descendants)
      !! Returns a logical mask over tips indicating descendants of a requested node.
      type(phylo_tree), intent(in) :: tree !! Rooted tree containing `node`.
      integer, intent(in) :: node !! Ape node number whose descendant tip set is requested.
      logical, allocatable, intent(out) :: descendants(:) !! Length-`n_tip` logical mask of descendant tips.
      integer, allocatable :: parent(:)
      integer :: current
      integer :: tip

      allocate(descendants(tree%n_tip))
      descendants = .false.
      if (.not. tree%valid()) return
      if (node < 1 .or. node > tree%total_nodes()) return
      if (node <= tree%n_tip) then
         descendants(node) = .true.
         return
      end if
      parent = parent_vector(tree)
      do tip = 1, tree%n_tip
         current = tip
         do
            if (current == node) then
               descendants(tip) = .true.
               exit
            end if
            if (parent(current) == 0) exit
            current = parent(current)
         end do
      end do
   end subroutine clade_tips


   pure subroutine topological_distance_ph85(tree_a, tree_b, distance, info)
      !! Computes ape's PH85 topological distance from canonical internal-node bipartitions.
      type(phylo_tree), intent(in) :: tree_a !! First tree; tip numbers must correspond to those in `tree_b`.
      type(phylo_tree), intent(in) :: tree_b !! Second tree with the same numbered tip set as `tree_a`.
      integer, intent(out) :: distance !! PH85 split distance; zero for identical canonical bipartition multisets.
      integer, intent(out) :: info !! Status code: zero on success, 1 for invalid trees, 2 for unequal tip counts.
      logical, allocatable :: splits_a(:, :)
      logical, allocatable :: splits_b(:, :)
      integer :: i
      integer :: j
      integer :: matches
      integer :: status

      distance = 0
      info = 0
      if (.not. tree_a%valid() .or. .not. tree_b%valid()) then
         info = 1
         return
      end if
      if (tree_a%n_tip /= tree_b%n_tip) then
         info = 2
         return
      end if
      call canonical_split_matrix(tree_a, splits_a, status)
      if (status /= 0) then
         info = 1
         return
      end if
      call canonical_split_matrix(tree_b, splits_b, status)
      if (status /= 0) then
         info = 1
         return
      end if
      matches = 0
      do i = 1, size(splits_b, 1)
         do j = 1, size(splits_a, 1)
            if (all(splits_b(i, :) .eqv. splits_a(j, :))) then
               matches = matches + 1
               exit
            end if
         end do
      end do
      distance = tree_a%n_node + tree_b%n_node - 2 * matches
   end subroutine topological_distance_ph85

   pure subroutine branch_score_distance(tree_a, tree_b, distance, info)
      !! Computes the Kuhner-Felsenstein branch-score distance used by `dist.topo(method="score")`.
      type(phylo_tree), intent(in) :: tree_a !! First tree with branch lengths and numbered tips matching `tree_b`.
      type(phylo_tree), intent(in) :: tree_b !! Second tree with branch lengths and the same numbered tip set.
      real(dp), intent(out) :: distance !! Square root of summed squared internal-split branch-length differences.
      integer, intent(out) :: info !! Status code: zero on success, 1 invalid tree, 2 unequal tips, 3 missing lengths.
      integer, allocatable :: edge_a(:)
      integer, allocatable :: edge_b(:)
      logical, allocatable :: splits_a(:, :)
      logical, allocatable :: splits_b(:, :)
      logical, allocatable :: matched_b(:)
      real(dp) :: sum_squares
      real(dp) :: length_a
      real(dp) :: length_b
      integer :: i
      integer :: j
      integer :: node_a
      integer :: node_b
      integer :: root_a
      integer :: root_b
      integer :: status
      logical :: found

      distance = 0.0_dp
      info = 0
      if (.not. tree_a%valid() .or. .not. tree_b%valid()) then
         info = 1
         return
      end if
      if (tree_a%n_tip /= tree_b%n_tip) then
         info = 2
         return
      end if
      if (.not. tree_a%has_lengths() .or. .not. tree_b%has_lengths()) then
         info = 3
         return
      end if
      call canonical_split_matrix(tree_a, splits_a, status)
      if (status /= 0) then
         info = 1
         return
      end if
      call canonical_split_matrix(tree_b, splits_b, status)
      if (status /= 0) then
         info = 1
         return
      end if
      edge_a = edge_index_to_child(tree_a)
      edge_b = edge_index_to_child(tree_b)
      allocate(matched_b(tree_b%n_node))
      matched_b = .false.
      root_a = tree_a%root()
      root_b = tree_b%root()
      sum_squares = 0.0_dp
      do i = 1, tree_a%n_node
         node_a = tree_a%n_tip + i
         if (node_a == root_a) cycle
         if (edge_a(node_a) == 0) then
            info = 1
            return
         end if
         length_a = tree_a%edge_length(edge_a(node_a))
         found = .false.
         do j = 1, tree_b%n_node
            node_b = tree_b%n_tip + j
            if (node_b == root_b .or. matched_b(j)) cycle
            if (.not. all(splits_a(i, :) .eqv. splits_b(j, :))) cycle
            if (edge_b(node_b) == 0) then
               info = 1
               return
            end if
            length_b = tree_b%edge_length(edge_b(node_b))
            sum_squares = sum_squares + (length_a - length_b)**2
            matched_b(j) = .true.
            found = .true.
            exit
         end do
         if (.not. found) sum_squares = sum_squares + length_a**2
      end do
      do j = 1, tree_b%n_node
         node_b = tree_b%n_tip + j
         if (node_b == root_b .or. matched_b(j)) cycle
         if (edge_b(node_b) == 0) then
            info = 1
            return
         end if
         length_b = tree_b%edge_length(edge_b(node_b))
         sum_squares = sum_squares + length_b**2
      end do
      distance = sqrt(sum_squares)
   end subroutine branch_score_distance

   pure subroutine canonical_split_matrix(tree, splits, info)
      !! Builds SHORTwise canonical internal-node bipartitions compatible with ape's PH85 comparison.
      type(phylo_tree), intent(in) :: tree !! Tree whose internal descendant-tip sets are canonicalized.
      logical, allocatable, intent(out) :: splits(:, :) !! Shape `(n_node,n_tip)` with shorter split sides selected.
      integer, intent(out) :: info !! Status code propagated from descendant-tip extraction.
      integer :: i
      integer :: n_selected

      call tip_descendant_matrix(tree, splits, info)
      if (info /= 0) return
      do i = 1, size(splits, 1)
         n_selected = count(splits(i, :))
         if (2 * n_selected > tree%n_tip) then
            splits(i, :) = .not. splits(i, :)
         else if (2 * n_selected == tree%n_tip) then
            if (.not. splits(i, 1)) splits(i, :) = .not. splits(i, :)
         end if
      end do
   end subroutine canonical_split_matrix

end module ape_topology
