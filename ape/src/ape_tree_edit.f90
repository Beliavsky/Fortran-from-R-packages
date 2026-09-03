! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
!
! Tree-editing algorithms translated from computational concepts in ape
! R/drop.tip.R, R/collapse.singles.R, R/multi2di.R, and R/root.R.
! Copyright holders and upstream provenance are documented in NOTICE.md.
module ape_tree_edit
   use r_kinds, only : dp
   use ape_types, only : phylo_tree, make_phylo_tree, parent_vector, child_counts, edge_index_to_child
   use ape_topology, only : clade_tips
   implicit none
   private

   public :: has_singles
   public :: is_rooted_tree
   public :: collapse_singles
   public :: drop_tips
   public :: keep_tips
   public :: extract_clade
   public :: reroot_node
   public :: root_outgroup
   public :: unroot_tree
   public :: di2multi
   public :: multi2di

contains

   pure logical function has_singles(tree) result(found)
      !! Tests whether any internal node has exactly one child, as in ape `has.singles`.
      type(phylo_tree), intent(in) :: tree !! Rooted tree whose internal out-degrees are inspected.
      integer, allocatable :: counts(:)
      integer :: node

      found = .false.
      if (.not. tree%valid()) return
      counts = child_counts(tree)
      do node = tree%n_tip + 1, tree%total_nodes()
         if (counts(node) == 1) then
            found = .true.
            return
         end if
      end do
   end function has_singles

   pure logical function is_rooted_tree(tree) result(rooted)
      !! Applies ape's structural rootedness test: the encoded root has at most two children.
      type(phylo_tree), intent(in) :: tree !! Tree whose root out-degree determines rooted versus unrooted encoding.
      integer, allocatable :: counts(:)
      integer :: root_node

      rooted = .false.
      if (.not. tree%valid()) return
      root_node = tree%root()
      if (root_node == 0) return
      counts = child_counts(tree)
      rooted = counts(root_node) <= 2
   end function is_rooted_tree

   pure subroutine collapse_singles(tree, collapsed, info)
      !! Removes internal nodes with one child and combines interior branch lengths along the collapsed path.
      type(phylo_tree), intent(in) :: tree !! Input rooted tree, with or without branch lengths.
      type(phylo_tree), intent(out) :: collapsed !! Tree after recursively removing every unary internal node.
      integer, intent(out) :: info !! Status code: zero on success, nonzero for an invalid or degenerate result.
      real(dp), allocatable :: length_to_parent(:)
      integer, allocatable :: parent(:)
      logical, allocatable :: active(:)

      info = 0
      if (.not. tree%valid()) then
         info = 1
         return
      end if
      call initialize_parent_state(tree, parent, length_to_parent, active)
      call collapse_active_singles(tree%n_tip, parent, length_to_parent, active)
      call rebuild_tree(tree%n_tip, parent, length_to_parent, active, collapsed, info, tree%has_lengths())
   end subroutine collapse_singles

   pure subroutine drop_tips(tree, tips, pruned, old_tip_numbers, info)
      !! Drops selected tip numbers, prunes empty ancestry, and collapses resulting singleton internal nodes.
      type(phylo_tree), intent(in) :: tree !! Input rooted tree whose tips use ape numbering `1:n_tip`.
      integer, intent(in) :: tips(:) !! Tip numbers to remove; duplicates are ignored and every value must be in range.
      type(phylo_tree), intent(out) :: pruned !! Pruned tree with kept tips renumbered in original ascending order.
      integer, allocatable, intent(out) :: old_tip_numbers(:) !! Original tip number corresponding to each new tip number.
      integer, intent(out) :: info !! Status code: zero on success, 1 invalid tree/tip index, 2 fewer than two tips remain.
      real(dp), allocatable :: length_to_parent(:)
      integer, allocatable :: parent(:)
      logical, allocatable :: active(:)
      logical, allocatable :: keep_tip(:)
      integer :: current
      integer :: i
      integer :: nkeep

      info = 0
      if (.not. tree%valid()) then
         info = 1
         return
      end if
      allocate(keep_tip(tree%n_tip))
      keep_tip = .true.
      do i = 1, size(tips)
         if (tips(i) < 1 .or. tips(i) > tree%n_tip) then
            info = 1
            return
         end if
         keep_tip(tips(i)) = .false.
      end do
      nkeep = count(keep_tip)
      if (nkeep < 2) then
         info = 2
         return
      end if
      allocate(old_tip_numbers(nkeep))
      old_tip_numbers = pack([(i, i = 1, tree%n_tip)], keep_tip)

      call initialize_parent_state(tree, parent, length_to_parent, active)
      active = .false.
      do i = 1, tree%n_tip
         if (.not. keep_tip(i)) cycle
         current = i
         do while (current /= 0)
            active(current) = .true.
            current = parent(current)
         end do
      end do
      call collapse_active_singles(tree%n_tip, parent, length_to_parent, active)
      call rebuild_tree(tree%n_tip, parent, length_to_parent, active, pruned, info, tree%has_lengths(), keep_tip)
   end subroutine drop_tips

   pure subroutine keep_tips(tree, tips, pruned, old_tip_numbers, info)
      !! Keeps selected tip numbers and drops all others, matching ape `keep.tip` semantics for numeric indices.
      type(phylo_tree), intent(in) :: tree !! Input rooted tree whose tips use ape numbering `1:n_tip`.
      integer, intent(in) :: tips(:) !! Unique or repeated tip numbers to retain; every value must be in range.
      type(phylo_tree), intent(out) :: pruned !! Tree induced by the requested tips after singleton collapse.
      integer, allocatable, intent(out) :: old_tip_numbers(:) !! Original tip number corresponding to each new tip number.
      integer, intent(out) :: info !! Status code propagated from validation or pruning.
      integer, allocatable :: drop(:)
      logical, allocatable :: keep(:)
      integer :: i
      integer :: ndrop

      info = 0
      if (.not. tree%valid()) then
         info = 1
         return
      end if
      allocate(keep(tree%n_tip))
      keep = .false.
      do i = 1, size(tips)
         if (tips(i) < 1 .or. tips(i) > tree%n_tip) then
            info = 1
            return
         end if
         keep(tips(i)) = .true.
      end do
      if (count(keep) < 2) then
         info = 2
         return
      end if
      ndrop = tree%n_tip - count(keep)
      allocate(drop(ndrop))
      drop = pack([(i, i = 1, tree%n_tip)], .not. keep)
      call drop_tips(tree, drop, pruned, old_tip_numbers, info)
   end subroutine keep_tips

   pure subroutine extract_clade(tree, node, clade, old_tip_numbers, info)
      !! Extracts the subtree induced by all tip descendants of an internal node.
      type(phylo_tree), intent(in) :: tree !! Rooted source tree containing the requested internal node.
      integer, intent(in) :: node !! Internal ape node number whose descendant clade is requested.
      type(phylo_tree), intent(out) :: clade !! Extracted tree with descendant tips renumbered from one.
      integer, allocatable, intent(out) :: old_tip_numbers(:) !! Original source-tree tip numbers in the extracted clade.
      integer, intent(out) :: info !! Status code: zero on success, 1 invalid node/tree, 2 fewer than two descendant tips.
      logical, allocatable :: descendants(:)
      integer, allocatable :: tips(:)
      integer :: i

      info = 0
      if (.not. tree%valid()) then
         info = 1
         return
      end if
      if (node <= tree%n_tip .or. node > tree%total_nodes()) then
         info = 1
         return
      end if
      call clade_tips(tree, node, descendants)
      if (count(descendants) < 2) then
         info = 2
         return
      end if
      allocate(tips(count(descendants)))
      tips = pack([(i, i = 1, tree%n_tip)], descendants)
      call keep_tips(tree, tips, clade, old_tip_numbers, info)
   end subroutine extract_clade

   pure subroutine reroot_node(tree, node, rooted_tree, info)
      !! Reroots a tree at an existing internal node by reversing edges on the old-root path.
      type(phylo_tree), intent(in) :: tree !! Input rooted tree; branch lengths, when present, are preserved on undirected edges.
      integer, intent(in) :: node !! Existing internal node that becomes the new root.
      type(phylo_tree), intent(out) :: rooted_tree !! Rerooted tree with the requested node renumbered as `n_tip + 1`.
      integer, intent(out) :: info !! Status code: zero on success, 1 invalid tree/node, 2 malformed root path.
      real(dp), allocatable :: length_to_parent(:)
      real(dp), allocatable :: new_length(:)
      integer, allocatable :: parent(:)
      integer, allocatable :: new_parent(:)
      logical, allocatable :: active(:)
      integer :: current
      integer :: next_node
      integer :: steps

      info = 0
      if (.not. tree%valid()) then
         info = 1
         return
      end if
      if (node <= tree%n_tip .or. node > tree%total_nodes()) then
         info = 1
         return
      end if
      call initialize_parent_state(tree, parent, length_to_parent, active)
      new_parent = parent
      new_length = length_to_parent
      current = node
      next_node = parent(current)
      new_parent(current) = 0
      steps = 0
      do while (next_node /= 0)
         steps = steps + 1
         if (steps > tree%total_nodes()) then
            info = 2
            return
         end if
         new_parent(next_node) = current
         new_length(next_node) = length_to_parent(current)
         current = next_node
         next_node = parent(current)
      end do
      call rebuild_tree(tree%n_tip, new_parent, new_length, active, rooted_tree, info, tree%has_lengths(), root_first=node)
   end subroutine reroot_node

   pure subroutine root_outgroup(tree, outgroup, rooted_tree, info, resolve_root)
      !! Reroots a tree with respect to a numeric outgroup, matching the computational semantics of ape `root.phylo`.
      type(phylo_tree), intent(in) :: tree !! Input tree whose tips use the same numeric identities referenced by `outgroup`.
      integer, intent(in) :: outgroup(:) !! Unique tips forming a monophyletic outgroup; nonempty and not all tips.
      type(phylo_tree), intent(out) :: rooted_tree !! Rerooted tree; root resolution may insert one zero-length basal edge.
      integer, intent(out) :: info !! Zero on success, 1 bad arguments, 2 nonmonophyly, or propagated edit status.
      logical, intent(in), optional :: resolve_root !! If true, add ape's zero branch below the ingroup MRCA to bifurcate the root.
      type(phylo_tree) :: base_tree
      type(phylo_tree) :: temporary
      logical, allocatable :: desired(:)
      logical, allocatable :: descendants(:)
      logical :: do_resolve
      integer :: e
      integer :: root_candidate
      integer :: unroot_info

      info = 0
      if (.not. tree%valid() .or. tree%n_tip < 2) then
         info = 1
         return
      end if
      if (.not. valid_tip_set(outgroup, tree%n_tip)) then
         info = 1
         return
      end if
      if (size(outgroup) >= tree%n_tip) then
         info = 1
         return
      end if
      do_resolve = .false.
      if (present(resolve_root)) do_resolve = resolve_root
      base_tree = tree
      if (size(outgroup) > 1 .and. is_rooted_tree(base_tree)) then
         call unroot_tree(base_tree, temporary, unroot_info)
         if (unroot_info /= 0) then
            info = unroot_info
            return
         end if
         base_tree = temporary
      end if

      allocate(desired(tree%n_tip))
      desired = .false.
      desired(outgroup) = .true.
      root_candidate = 0
      do e = 1, base_tree%nedge()
         call clade_tips(base_tree, base_tree%edge(e, 2), descendants)
         if (all(descendants .eqv. desired)) then
            if (base_tree%edge(e, 1) > base_tree%n_tip) then
               root_candidate = base_tree%edge(e, 1)
            else
               root_candidate = base_tree%edge(e, 2)
            end if
            exit
         end if
         if (all((.not. descendants) .eqv. desired)) then
            if (base_tree%edge(e, 2) > base_tree%n_tip) then
               root_candidate = base_tree%edge(e, 2)
            else
               root_candidate = base_tree%edge(e, 1)
            end if
            exit
         end if
      end do
      if (root_candidate <= base_tree%n_tip) then
         info = 2
         return
      end if
      call reroot_node(base_tree, root_candidate, rooted_tree, info)
      if (info /= 0 .or. .not. do_resolve) return
      call resolve_outgroup_root(rooted_tree, desired, temporary, info)
      if (info == 0) rooted_tree = temporary
   end subroutine root_outgroup

   pure subroutine unroot_tree(tree, unrooted, info)
      !! Converts a bifurcating rooted encoding into ape's degree-three-root unrooted representation.
      type(phylo_tree), intent(in) :: tree !! Rooted tree with at least three tips and a root having exactly two children.
      type(phylo_tree), intent(out) :: unrooted !! Unrooted representation with one fewer internal node.
      integer, intent(out) :: info !! Status code: zero on success; nonzero for invalid, already-unrooted, or degenerate trees.
      real(dp), allocatable :: length_to_parent(:)
      integer, allocatable :: counts(:)
      integer, allocatable :: parent(:)
      logical, allocatable :: active(:)
      integer :: child_a
      integer :: child_b
      integer :: e
      integer :: new_root
      integer :: other
      integer :: root_node

      info = 0
      if (.not. tree%valid() .or. tree%n_tip < 3) then
         info = 1
         return
      end if
      root_node = tree%root()
      counts = child_counts(tree)
      if (counts(root_node) /= 2) then
         info = 2
         return
      end if
      child_a = 0
      child_b = 0
      do e = 1, tree%nedge()
         if (tree%edge(e, 1) /= root_node) cycle
         if (child_a == 0) then
            child_a = tree%edge(e, 2)
         else
            child_b = tree%edge(e, 2)
         end if
      end do
      if (child_a == 0 .or. child_b == 0) then
         info = 1
         return
      end if
      if (child_a > tree%n_tip) then
         new_root = child_a
         other = child_b
      else if (child_b > tree%n_tip) then
         new_root = child_b
         other = child_a
      else
         info = 3
         return
      end if

      call initialize_parent_state(tree, parent, length_to_parent, active)
      parent(new_root) = 0
      parent(other) = new_root
      if (tree%has_lengths()) length_to_parent(other) = length_to_parent(child_a) + length_to_parent(child_b)
      active(root_node) = .false.
      call rebuild_tree(tree%n_tip, parent, length_to_parent, active, unrooted, info, tree%has_lengths(), root_first=new_root)
   end subroutine unroot_tree

   pure subroutine di2multi(tree, tol, multifurcating, info)
      !! Contracts internal edges shorter than `tol`, translating ape `di2multi(..., tip2root=FALSE)`.
      type(phylo_tree), intent(in) :: tree !! Input tree with branch lengths.
      real(dp), intent(in) :: tol !! Nonnegative threshold; internal edges with length strictly below it are contracted.
      type(phylo_tree), intent(out) :: multifurcating !! Tree after all qualifying internal-edge contractions.
      integer, intent(out) :: info !! Status code: zero on success, nonzero for invalid input or negative tolerance.
      real(dp), allocatable :: length_to_parent(:)
      integer, allocatable :: parent(:)
      logical, allocatable :: active(:)
      integer :: node
      integer :: p
      logical :: changed

      info = 0
      if (.not. tree%valid() .or. .not. tree%has_lengths() .or. tol < 0.0_dp) then
         info = 1
         return
      end if
      call initialize_parent_state(tree, parent, length_to_parent, active)
      do
         changed = .false.
         do node = tree%n_tip + 1, tree%total_nodes()
            if (.not. active(node)) cycle
            p = parent(node)
            if (p == 0) cycle
            if (length_to_parent(node) >= tol) cycle
            call contract_node(node, parent, active)
            changed = .true.
         end do
         if (.not. changed) exit
      end do
      call rebuild_tree(tree%n_tip, parent, length_to_parent, active, multifurcating, info, .true.)
   end subroutine di2multi

   pure subroutine multi2di(tree, dichotomous, info)
      !! Deterministically resolves multifurcations with zero-length internal branches, matching `multi2di(random=FALSE)`.
      type(phylo_tree), intent(in) :: tree !! Input rooted or unrooted tree, with or without branch lengths.
      type(phylo_tree), intent(out) :: dichotomous !! Tree in which every internal node has at most two children.
      integer, intent(out) :: info !! Status code: zero on success, nonzero for an invalid tree.
      integer, allocatable :: counts(:)
      integer, allocatable :: child(:)
      integer, allocatable :: child_edge(:)
      integer, allocatable :: edges(:, :)
      real(dp), allocatable :: lengths(:)
      integer :: added
      integer :: e
      integer :: edge_pos
      integer :: i
      integer :: j
      integer :: n_children
      integer :: new_node
      integer :: next_node
      integer :: node
      integer :: old_total
      integer :: parent_node
      logical :: with_lengths

      info = 0
      if (.not. tree%valid()) then
         info = 1
         return
      end if
      counts = child_counts(tree)
      added = 0
      do node = tree%n_tip + 1, tree%total_nodes()
         if (counts(node) > 2) added = added + counts(node) - 2
      end do
      if (added == 0) then
         dichotomous = tree
         return
      end if
      allocate(edges(tree%nedge() + added, 2))
      allocate(lengths(tree%nedge() + added))
      lengths = 0.0_dp
      edge_pos = 0
      old_total = tree%total_nodes()
      next_node = old_total + 1
      with_lengths = tree%has_lengths()

      do node = 1, tree%total_nodes()
         n_children = counts(node)
         if (n_children == 0) cycle
         allocate(child(n_children), child_edge(n_children))
         j = 0
         do e = 1, tree%nedge()
            if (tree%edge(e, 1) /= node) cycle
            j = j + 1
            child(j) = tree%edge(e, 2)
            child_edge(j) = e
         end do
         if (n_children <= 2) then
            do i = 1, n_children
               edge_pos = edge_pos + 1
               edges(edge_pos, :) = [node, child(i)]
               if (with_lengths) lengths(edge_pos) = tree%edge_length(child_edge(i))
            end do
         else
            parent_node = node
            do i = 1, n_children - 2
               edge_pos = edge_pos + 1
               edges(edge_pos, :) = [parent_node, child(i)]
               if (with_lengths) lengths(edge_pos) = tree%edge_length(child_edge(i))
               new_node = next_node
               next_node = next_node + 1
               edge_pos = edge_pos + 1
               edges(edge_pos, :) = [parent_node, new_node]
               lengths(edge_pos) = 0.0_dp
               parent_node = new_node
            end do
            edge_pos = edge_pos + 1
            edges(edge_pos, :) = [parent_node, child(n_children - 1)]
            if (with_lengths) lengths(edge_pos) = tree%edge_length(child_edge(n_children - 1))
            edge_pos = edge_pos + 1
            edges(edge_pos, :) = [parent_node, child(n_children)]
            if (with_lengths) lengths(edge_pos) = tree%edge_length(child_edge(n_children))
         end if
         deallocate(child, child_edge)
      end do
      if (with_lengths) then
         call renumber_tree(tree%n_tip, edges, lengths, dichotomous, info)
      else
         call renumber_tree(tree%n_tip, edges, result=dichotomous, info=info)
      end if
   end subroutine multi2di

   pure subroutine initialize_parent_state(tree, parent, length_to_parent, active)
      !! Converts an edge-list tree into mutable parent and incoming-length vectors.
      type(phylo_tree), intent(in) :: tree !! Valid source tree.
      integer, allocatable, intent(out) :: parent(:) !! Immediate parent of every source node, zero at the root.
      real(dp), allocatable, intent(out) :: length_to_parent(:) !! Incoming branch length for each node, zero when unavailable/root.
      logical, allocatable, intent(out) :: active(:) !! Activity mask initialized true for every source node.
      integer, allocatable :: edge_for_child(:)
      integer :: node

      parent = parent_vector(tree)
      allocate(length_to_parent(tree%total_nodes()), active(tree%total_nodes()))
      length_to_parent = 0.0_dp
      active = .true.
      if (.not. tree%has_lengths()) return
      edge_for_child = edge_index_to_child(tree)
      do node = 1, tree%total_nodes()
         if (edge_for_child(node) /= 0) length_to_parent(node) = tree%edge_length(edge_for_child(node))
      end do
   end subroutine initialize_parent_state

   pure subroutine collapse_active_singles(n_tip, parent, length_to_parent, active)
      !! Recursively removes unary active internal nodes from a mutable parent representation.
      integer, intent(in) :: n_tip !! Number of original tips; nodes above this value are internal.
      integer, intent(inout) :: parent(:) !! Mutable immediate-parent vector.
      real(dp), intent(inout) :: length_to_parent(:) !! Mutable incoming branch lengths, summed across interior collapses.
      logical, intent(inout) :: active(:) !! Mutable node activity mask.
      integer :: child
      integer :: child_count
      integer :: node
      integer :: p
      logical :: changed

      do
         changed = .false.
         do node = n_tip + 1, size(active)
            if (.not. active(node)) cycle
            child = 0
            child_count = 0
            call active_child_summary(node, parent, active, child_count, child)
            if (child_count /= 1) cycle
            p = parent(node)
            if (p == 0) then
               parent(child) = 0
               length_to_parent(child) = 0.0_dp
            else
               parent(child) = p
               length_to_parent(child) = length_to_parent(child) + length_to_parent(node)
            end if
            active(node) = .false.
            parent(node) = 0
            length_to_parent(node) = 0.0_dp
            changed = .true.
         end do
         if (.not. changed) exit
      end do
   end subroutine collapse_active_singles

   pure subroutine active_child_summary(node, parent, active, child_count, only_child)
      !! Counts active children of a node and returns the child when the count is one.
      integer, intent(in) :: node !! Candidate parent node.
      integer, intent(in) :: parent(:) !! Immediate-parent vector.
      logical, intent(in) :: active(:) !! Node activity mask.
      integer, intent(out) :: child_count !! Number of active nodes whose immediate parent is `node`.
      integer, intent(out) :: only_child !! Sole child when `child_count=1`, otherwise the last child encountered or zero.
      integer :: i

      child_count = 0
      only_child = 0
      do i = 1, size(parent)
         if (.not. active(i)) cycle
         if (parent(i) /= node) cycle
         child_count = child_count + 1
         only_child = i
      end do
   end subroutine active_child_summary

   pure subroutine contract_node(node, parent, active)
      !! Contracts one active internal node by attaching all of its children directly to its parent.
      integer, intent(in) :: node !! Nonroot internal node to remove.
      integer, intent(inout) :: parent(:) !! Mutable immediate-parent vector.
      logical, intent(inout) :: active(:) !! Mutable activity mask; `node` is deactivated.
      integer :: i
      integer :: p

      p = parent(node)
      do i = 1, size(parent)
         if (.not. active(i)) cycle
         if (parent(i) == node) parent(i) = p
      end do
      active(node) = .false.
      parent(node) = 0
   end subroutine contract_node

   pure subroutine rebuild_tree(old_n_tip, parent, length_to_parent, active, result, info, with_lengths, keep_tip, root_first)
      !! Rebuilds and renumbers a tree from active parent relationships.
      integer, intent(in) :: old_n_tip !! Tip-number boundary in the mutable source representation.
      integer, intent(in) :: parent(:) !! Parent vector over source node numbers.
      real(dp), intent(in) :: length_to_parent(:) !! Incoming branch lengths over source node numbers.
      logical, intent(in) :: active(:) !! Mask of source nodes retained in the rebuilt tree.
      type(phylo_tree), intent(out) :: result !! Rebuilt tree with contiguous ape-style tip/internal numbering.
      integer, intent(out) :: info !! Status code: zero on success, nonzero if the active relationships are degenerate.
      logical, intent(in) :: with_lengths !! Whether to attach branch lengths to the result.
      logical, intent(in), optional :: keep_tip(:) !! Optional old-tip mask; absent means every active old tip is retained.
      integer, intent(in), optional :: root_first !! Optional active internal node forced to become new node `n_tip + 1`.
      integer, allocatable :: edges(:, :)
      real(dp), allocatable :: lengths(:)
      integer, allocatable :: map(:)
      logical, allocatable :: tip_mask(:)
      integer :: e
      integer :: i
      integer :: n_edge
      integer :: n_internal
      integer :: n_tip
      integer :: next_internal
      integer :: p
      integer :: root_node

      info = 0
      allocate(tip_mask(old_n_tip))
      if (present(keep_tip)) then
         if (size(keep_tip) /= old_n_tip) then
            info = 1
            return
         end if
         tip_mask = keep_tip .and. active(1:old_n_tip)
      else
         tip_mask = active(1:old_n_tip)
      end if
      n_tip = count(tip_mask)
      if (n_tip < 2) then
         info = 2
         return
      end if
      root_node = 0
      do i = old_n_tip + 1, size(active)
         if (.not. active(i)) cycle
         if (parent(i) == 0) then
            if (root_node /= 0) then
               info = 3
               return
            end if
            root_node = i
         end if
      end do
      if (present(root_first)) then
         if (root_first < old_n_tip + 1 .or. root_first > size(active)) then
            info = 3
            return
         end if
         if (.not. active(root_first) .or. parent(root_first) /= 0) then
            info = 3
            return
         end if
         root_node = root_first
      end if
      if (root_node == 0) then
         info = 3
         return
      end if

      allocate(map(size(active)))
      map = 0
      e = 0
      do i = 1, old_n_tip
         if (.not. tip_mask(i)) cycle
         e = e + 1
         map(i) = e
      end do
      n_internal = count(active(old_n_tip + 1:))
      if (n_internal < 1) then
         info = 3
         return
      end if
      map(root_node) = n_tip + 1
      next_internal = n_tip + 2
      do i = old_n_tip + 1, size(active)
         if (.not. active(i) .or. i == root_node) cycle
         map(i) = next_internal
         next_internal = next_internal + 1
      end do

      n_edge = count(active) - 1
      if (n_edge /= n_tip + n_internal - 1) then
         info = 3
         return
      end if
      allocate(edges(n_edge, 2), lengths(n_edge))
      lengths = 0.0_dp
      e = 0
      do i = 1, size(active)
         if (.not. active(i) .or. i == root_node) cycle
         p = parent(i)
         if (p <= 0 .or. p > size(active) .or. .not. active(p)) then
            info = 3
            return
         end if
         if (map(i) == 0 .or. map(p) == 0) then
            info = 3
            return
         end if
         e = e + 1
         edges(e, :) = [map(p), map(i)]
         if (with_lengths) lengths(e) = length_to_parent(i)
      end do
      if (with_lengths) then
         result = make_phylo_tree(n_tip, edges, lengths)
      else
         result = make_phylo_tree(n_tip, edges)
      end if
      if (.not. result%valid()) info = 4
   end subroutine rebuild_tree

   pure logical function valid_tip_set(tips, n_tip) result(ok)
      !! Checks that a numeric tip set is nonempty, in range, and contains no duplicates.
      integer, intent(in) :: tips(:) !! Candidate tip indices.
      integer, intent(in) :: n_tip !! Number of available tips, defining the valid range `1:n_tip`.
      logical, allocatable :: seen(:)
      integer :: i

      ok = .false.
      if (size(tips) == 0 .or. n_tip < 1) return
      allocate(seen(n_tip))
      seen = .false.
      do i = 1, size(tips)
         if (tips(i) < 1 .or. tips(i) > n_tip) return
         if (seen(tips(i))) return
         seen(tips(i)) = .true.
      end do
      ok = .true.
   end function valid_tip_set

   pure subroutine resolve_outgroup_root(tree, outgroup_mask, resolved, info)
      !! Resolves a basal multifurcation by grouping all ingroup root branches below one new zero-length internal edge.
      type(phylo_tree), intent(in) :: tree !! Rerooted tree whose root is adjacent to the monophyletic outgroup branch.
      logical, intent(in) :: outgroup_mask(:) !! Tip-membership mask of the requested outgroup in current tip numbering.
      type(phylo_tree), intent(out) :: resolved !! Root-bifurcating tree with unchanged tip-to-tip path lengths.
      integer, intent(out) :: info !! Zero on success, 1 if no outgroup branch is found, otherwise renumbering status.
      logical, allocatable :: descendants(:)
      integer, allocatable :: edge(:, :)
      real(dp), allocatable :: edge_length(:)
      integer :: e
      integer :: outgroup_edge
      integer :: root_node
      integer :: new_node

      info = 0
      if (.not. tree%valid()) then
         info = 1
         return
      end if
      root_node = tree%root()
      if (count(tree%edge(:, 1) == root_node) <= 2) then
         resolved = tree
         return
      end if
      outgroup_edge = 0
      do e = 1, tree%nedge()
         if (tree%edge(e, 1) /= root_node) cycle
         call clade_tips(tree, tree%edge(e, 2), descendants)
         if (all(descendants .eqv. outgroup_mask)) then
            outgroup_edge = e
            exit
         end if
      end do
      if (outgroup_edge == 0) then
         info = 1
         return
      end if
      new_node = tree%total_nodes() + 1
      allocate(edge(tree%nedge() + 1, 2))
      edge(1:tree%nedge(), :) = tree%edge
      do e = 1, tree%nedge()
         if (edge(e, 1) == root_node .and. e /= outgroup_edge) edge(e, 1) = new_node
      end do
      edge(tree%nedge() + 1, :) = [root_node, new_node]
      if (tree%has_lengths()) then
         allocate(edge_length(tree%nedge() + 1))
         edge_length(1:tree%nedge()) = tree%edge_length
         edge_length(tree%nedge() + 1) = 0.0_dp
         call renumber_tree(tree%n_tip, edge, edge_length, resolved, info)
      else
         call renumber_tree(tree%n_tip, edge, result=resolved, info=info)
      end if
   end subroutine resolve_outgroup_root

   pure subroutine renumber_tree(n_tip, edge, edge_length, result, info)
      !! Renumbers possibly sparse internal node identifiers into a contiguous ape-style range.
      integer, intent(in) :: n_tip !! Number of terminal taxa, already numbered `1:n_tip` in `edge`.
      integer, intent(in) :: edge(:, :) !! Directed edge matrix with arbitrary positive internal node numbers.
      real(dp), intent(in), optional :: edge_length(:) !! Optional branch lengths aligned with edge rows.
      type(phylo_tree), intent(out) :: result !! Renumbered tree with the root assigned `n_tip + 1`.
      integer, intent(out) :: info !! Status code: zero on success, nonzero for malformed topology or lengths.
      integer, allocatable :: map(:)
      integer, allocatable :: new_edge(:, :)
      logical, allocatable :: is_child(:)
      logical, allocatable :: is_internal(:)
      integer :: i
      integer :: max_node
      integer :: next_node
      integer :: root_node

      info = 0
      if (size(edge, 2) /= 2 .or. size(edge, 1) == 0) then
         info = 1
         return
      end if
      if (present(edge_length)) then
         if (size(edge_length) /= size(edge, 1)) then
            info = 1
            return
         end if
      end if
      max_node = maxval(edge)
      allocate(map(max_node), is_child(max_node), is_internal(max_node))
      map = 0
      is_child = .false.
      is_internal = .false.
      do i = 1, n_tip
         map(i) = i
      end do
      do i = 1, size(edge, 1)
         if (edge(i, 1) > n_tip) is_internal(edge(i, 1)) = .true.
         if (edge(i, 2) > n_tip) is_internal(edge(i, 2)) = .true.
         is_child(edge(i, 2)) = .true.
      end do
      root_node = 0
      do i = n_tip + 1, max_node
         if (.not. is_internal(i) .or. is_child(i)) cycle
         if (root_node /= 0) then
            info = 2
            return
         end if
         root_node = i
      end do
      if (root_node == 0) then
         info = 2
         return
      end if
      map(root_node) = n_tip + 1
      next_node = n_tip + 2
      do i = n_tip + 1, max_node
         if (.not. is_internal(i) .or. i == root_node) cycle
         map(i) = next_node
         next_node = next_node + 1
      end do
      allocate(new_edge(size(edge, 1), 2))
      do i = 1, size(edge, 1)
         if (edge(i, 1) > max_node .or. edge(i, 2) > max_node) then
            info = 2
            return
         end if
         if (map(edge(i, 1)) == 0 .or. map(edge(i, 2)) == 0) then
            info = 2
            return
         end if
         new_edge(i, :) = [map(edge(i, 1)), map(edge(i, 2))]
      end do
      if (present(edge_length)) then
         result = make_phylo_tree(n_tip, new_edge, edge_length)
      else
         result = make_phylo_tree(n_tip, new_edge)
      end if
      if (.not. result%valid()) info = 3
   end subroutine renumber_tree

end module ape_tree_edit
