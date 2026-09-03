! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
!
! Modern Fortran data structures derived from computational concepts in
! R package ape 5.8-1. This translation is not part of upstream ape.
module ape_types
   use r_kinds, only : dp
   implicit none
   private

   type, public :: phylo_tree
      integer :: n_tip = 0
      integer :: n_node = 0
      integer, allocatable :: edge(:, :)
      real(dp), allocatable :: edge_length(:)
   contains
      procedure :: nedge => tree_nedge
      procedure :: total_nodes => tree_total_nodes
      procedure :: root => tree_root
      procedure :: has_lengths => tree_has_lengths
      procedure :: valid => tree_valid
   end type phylo_tree

   public :: make_phylo_tree
   public :: parent_vector
   public :: child_counts
   public :: edge_index_to_child

contains

   pure function make_phylo_tree(n_tip, edge, edge_length) result(tree)
      !! Constructs a phylogenetic tree from an ape-style edge matrix.
      integer, intent(in) :: n_tip !! Number of terminal taxa; must be positive.
      integer, intent(in) :: edge(:, :) !! Edge matrix with shape `(nedge, 2)` using 1-based node numbers.
      real(dp), intent(in), optional :: edge_length(:) !! Optional branch lengths with one value per edge.
      type(phylo_tree) :: tree
      integer :: max_node

      tree%n_tip = n_tip
      allocate(tree%edge(size(edge, 1), 2))
      tree%edge = edge
      if (size(edge, 1) == 0) then
         tree%n_node = 0
      else
         max_node = maxval(edge)
         tree%n_node = max(0, max_node - n_tip)
      end if
      if (present(edge_length)) then
         allocate(tree%edge_length(size(edge_length)))
         tree%edge_length = edge_length
      end if
   end function make_phylo_tree

   pure elemental integer function tree_nedge(self) result(value)
      !! Returns the number of directed edges stored in the tree.
      class(phylo_tree), intent(in) :: self !! Tree whose edge count is requested.

      if (allocated(self%edge)) then
         value = size(self%edge, 1)
      else
         value = 0
      end if
   end function tree_nedge

   pure elemental integer function tree_total_nodes(self) result(value)
      !! Returns the total number of tip and internal nodes.
      class(phylo_tree), intent(in) :: self !! Tree whose node count is requested.

      value = self%n_tip + self%n_node
   end function tree_total_nodes

   pure elemental integer function tree_root(self) result(root_node)
      !! Returns the unique root node number, or zero for an invalid/unrooted encoding.
      class(phylo_tree), intent(in) :: self !! Tree encoded by parent-child directed edges.
      logical, allocatable :: is_child(:)
      integer :: i
      integer :: n_total
      integer :: candidates

      root_node = 0
      n_total = self%total_nodes()
      if (.not. allocated(self%edge) .or. n_total <= 0) return
      allocate(is_child(n_total))
      is_child = .false.
      do i = 1, size(self%edge, 1)
         if (self%edge(i, 2) >= 1 .and. self%edge(i, 2) <= n_total) then
            is_child(self%edge(i, 2)) = .true.
         end if
      end do
      candidates = 0
      do i = self%n_tip + 1, n_total
         if (.not. is_child(i)) then
            root_node = i
            candidates = candidates + 1
         end if
      end do
      if (candidates /= 1) root_node = 0
   end function tree_root

   pure elemental logical function tree_has_lengths(self) result(has_lengths)
      !! Reports whether a branch length is present for every edge.
      class(phylo_tree), intent(in) :: self !! Tree whose branch-length storage is checked.

      has_lengths = allocated(self%edge_length)
      if (has_lengths) has_lengths = size(self%edge_length) == self%nedge()
   end function tree_has_lengths

   pure elemental logical function tree_valid(self) result(ok)
      !! Checks basic ape-style tree invariants without requiring a binary topology.
      class(phylo_tree), intent(in) :: self !! Tree whose shape, node numbering, and parent-child encoding are validated.
      integer :: i
      integer :: n_total
      integer :: root_node
      integer, allocatable :: parents(:)

      ok = .false.
      if (self%n_tip < 1 .or. self%n_node < 0) return
      if (.not. allocated(self%edge)) return
      if (size(self%edge, 2) /= 2) return
      if (self%has_lengths()) then
         if (size(self%edge_length) /= size(self%edge, 1)) return
      else if (allocated(self%edge_length)) then
         return
      end if
      n_total = self%total_nodes()
      if (size(self%edge, 1) > 0) then
         if (minval(self%edge) < 1 .or. maxval(self%edge) > n_total) return
      end if
      allocate(parents(n_total))
      parents = 0
      do i = 1, size(self%edge, 1)
         if (self%edge(i, 1) == self%edge(i, 2)) return
         if (parents(self%edge(i, 2)) /= 0) return
         parents(self%edge(i, 2)) = self%edge(i, 1)
      end do
      root_node = self%root()
      if (root_node == 0) return
      if (parents(root_node) /= 0) return
      do i = 1, n_total
         if (i == root_node) cycle
         if (parents(i) == 0) return
      end do
      ok = .true.
   end function tree_valid

   pure function parent_vector(tree) result(parent)
      !! Builds a node-indexed lookup giving each node's immediate parent.
      type(phylo_tree), intent(in) :: tree !! Tree whose directed edge matrix defines the parent relationships.
      integer, allocatable :: parent(:)
      integer :: i
      integer :: n_total

      n_total = tree%total_nodes()
      allocate(parent(n_total))
      parent = 0
      if (.not. allocated(tree%edge)) return
      do i = 1, size(tree%edge, 1)
         if (tree%edge(i, 2) >= 1 .and. tree%edge(i, 2) <= n_total) then
            parent(tree%edge(i, 2)) = tree%edge(i, 1)
         end if
      end do
   end function parent_vector

   pure function child_counts(tree) result(counts)
      !! Counts outgoing child edges for every tip and internal node.
      type(phylo_tree), intent(in) :: tree !! Tree whose node out-degrees are requested.
      integer, allocatable :: counts(:)
      integer :: i

      allocate(counts(tree%total_nodes()))
      counts = 0
      if (.not. allocated(tree%edge)) return
      do i = 1, size(tree%edge, 1)
         counts(tree%edge(i, 1)) = counts(tree%edge(i, 1)) + 1
      end do
   end function child_counts

   pure function edge_index_to_child(tree) result(index_by_child)
      !! Maps each non-root node to the edge index on which it occurs as a descendant.
      type(phylo_tree), intent(in) :: tree !! Tree whose edge lookup is constructed.
      integer, allocatable :: index_by_child(:)
      integer :: i

      allocate(index_by_child(tree%total_nodes()))
      index_by_child = 0
      if (.not. allocated(tree%edge)) return
      do i = 1, size(tree%edge, 1)
         index_by_child(tree%edge(i, 2)) = i
      end do
   end function edge_index_to_child

end module ape_types
