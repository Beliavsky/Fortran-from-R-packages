! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
!
! Split and bipartition collection algorithms translated from ape
! R/dist.topo.R and src/bitsplits.c. Copyright holders and detailed
! upstream provenance are documented in NOTICE.md.
module ape_splits
   use r_kinds, only : dp
   use ape_types, only : phylo_tree, make_phylo_tree
   use ape_topology, only : clade_tips
   implicit none
   private

   type, public :: split_collection
      integer :: n_tip = 0
      integer :: n_split = 0
      logical, allocatable :: split(:, :)
      integer, allocatable :: frequency(:)
   end type split_collection

   public :: prop_part
   public :: bitsplits
   public :: count_bipartitions
   public :: consensus_tree
   public :: prop_clades
   public :: tree_bipartitions

contains

   pure subroutine prop_part(trees, collection, info)
      !! Collects rooted descendant-tip clades and their frequencies, analogous to ape `prop.part`.
      type(phylo_tree), intent(in) :: trees(:) !! Trees with common tip numbering; each internal clade contributes once per tree.
      type(split_collection), intent(out) :: collection !! Unique rooted clades and occurrence counts across the input trees.
      integer, intent(out) :: info !! Status code: zero on success, 1 for invalid trees, 2 for inconsistent tip counts.
      logical, allocatable :: mask(:)
      logical, allocatable :: storage(:, :)
      integer, allocatable :: frequency(:)
      integer :: i
      integer :: max_split
      integer :: node
      integer :: n_tip
      integer :: used

      collection%n_tip = 0
      collection%n_split = 0
      info = 0
      if (size(trees) == 0) then
         allocate(collection%split(0, 0), collection%frequency(0))
         return
      end if
      if (.not. trees(1)%valid()) then
         info = 1
         return
      end if
      n_tip = trees(1)%n_tip
      max_split = 0
      do i = 1, size(trees)
         if (.not. trees(i)%valid()) then
            info = 1
            return
         end if
         if (trees(i)%n_tip /= n_tip) then
            info = 2
            return
         end if
         max_split = max_split + trees(i)%n_node
      end do
      allocate(storage(max_split, n_tip), frequency(max_split))
      storage = .false.
      frequency = 0
      used = 0
      do i = 1, size(trees)
         do node = trees(i)%n_tip + 1, trees(i)%total_nodes()
            call clade_tips(trees(i), node, mask)
            call store_split(mask, storage, frequency, used)
         end do
      end do
      call finish_collection(n_tip, storage, frequency, used, collection)
   end subroutine prop_part

   pure subroutine bitsplits(trees, collection, info)
      !! Collects unique nontrivial unrooted bipartitions with ape `bitsplits` one-wise canonicalization.
      type(phylo_tree), intent(in) :: trees(:) !! Trees with common tip numbering; rooted storage is treated as undirected.
      type(split_collection), intent(out) :: collection !! Unique nontrivial splits and their frequencies across trees.
      integer, intent(out) :: info !! Status code: zero on success, 1 for invalid trees, 2 for inconsistent tip counts.
      logical, allocatable :: one_tree(:, :)
      logical, allocatable :: storage(:, :)
      integer, allocatable :: frequency(:)
      integer :: i
      integer :: j
      integer :: max_split
      integer :: n_tip
      integer :: split_info
      integer :: used

      collection%n_tip = 0
      collection%n_split = 0
      info = 0
      if (size(trees) == 0) then
         allocate(collection%split(0, 0), collection%frequency(0))
         return
      end if
      if (.not. trees(1)%valid()) then
         info = 1
         return
      end if
      n_tip = trees(1)%n_tip
      max_split = 0
      do i = 1, size(trees)
         if (.not. trees(i)%valid()) then
            info = 1
            return
         end if
         if (trees(i)%n_tip /= n_tip) then
            info = 2
            return
         end if
         max_split = max_split + max(0, trees(i)%nedge() - n_tip)
      end do
      allocate(storage(max_split, n_tip), frequency(max_split))
      storage = .false.
      frequency = 0
      used = 0
      do i = 1, size(trees)
         call tree_bipartitions(trees(i), one_tree, split_info)
         if (split_info /= 0) then
            info = split_info
            return
         end if
         do j = 1, size(one_tree, 1)
            call store_split(one_tree(j, :), storage, frequency, used)
         end do
      end do
      call finish_collection(n_tip, storage, frequency, used, collection)
   end subroutine bitsplits

   pure subroutine count_bipartitions(tree, reference_trees, counts, info)
      !! Counts how often each nontrivial split of one tree occurs among reference trees, analogous to ape `countBipartitions`.
      type(phylo_tree), intent(in) :: tree !! Query tree whose nontrivial bipartitions are counted.
      type(phylo_tree), intent(in) :: reference_trees(:) !! Reference trees with the same tip numbering as the query tree.
      integer, allocatable, intent(out) :: counts(:) !! Frequency of each query split in the order returned by `tree_bipartitions`.
      integer, intent(out) :: info !! Status code propagated from split extraction and collection.
      type(split_collection) :: reference
      logical, allocatable :: query(:, :)
      integer :: i
      integer :: j

      call tree_bipartitions(tree, query, info)
      if (info /= 0) return
      call bitsplits(reference_trees, reference, info)
      if (info /= 0) return
      allocate(counts(size(query, 1)))
      counts = 0
      do i = 1, size(query, 1)
         do j = 1, reference%n_split
            if (all(query(i, :) .eqv. reference%split(j, :))) then
               counts(i) = reference%frequency(j)
               exit
            end if
         end do
      end do
   end subroutine count_bipartitions

   pure subroutine prop_clades(tree, reference_trees, counts, info, rooted)
      !! Counts support for each internal-node clade, matching ape `prop.clades` rooted or SHORTwise unrooted semantics.
      type(phylo_tree), intent(in) :: tree !! Query tree whose internal-node clades are assigned support counts.
      type(phylo_tree), intent(in) :: reference_trees(:) !! Reference trees sharing the query tip numbering.
      integer, allocatable, intent(out) :: counts(:) !! Internal-node support counts; -1 marks clades absent from references.
      integer, intent(out) :: info !! Status code: zero on success, 1 for invalid trees, or 2 for inconsistent tip counts.
      logical, intent(in), optional :: rooted !! True uses rooted clades; false uses ape SHORTwise split canonicalization.
      type(split_collection) :: reference
      logical, allocatable :: mask(:)
      integer :: i
      integer :: j
      integer :: node
      integer :: n_tip
      logical :: rooted_mode

      info = 0
      if (.not. tree%valid()) then
         info = 1
         return
      end if
      n_tip = tree%n_tip
      rooted_mode = .false.
      if (present(rooted)) rooted_mode = rooted
      if (size(reference_trees) == 0) then
         allocate(counts(tree%n_node))
         counts = -1
         return
      end if
      do i = 1, size(reference_trees)
         if (.not. reference_trees(i)%valid()) then
            info = 1
            return
         end if
         if (reference_trees(i)%n_tip /= n_tip) then
            info = 2
            return
         end if
      end do
      if (rooted_mode) then
         call prop_part(reference_trees, reference, info)
      else
         call shortwise_part(reference_trees, reference, info)
      end if
      if (info /= 0) return
      allocate(counts(tree%n_node))
      counts = -1
      do i = 1, tree%n_node
         node = n_tip + i
         call clade_tips(tree, node, mask)
         if (.not. rooted_mode) call shortwise_mask(mask)
         do j = 1, reference%n_split
            if (all(mask .eqv. reference%split(j, :))) then
               counts(i) = reference%frequency(j)
               exit
            end if
         end do
      end do
   end subroutine prop_clades

   pure subroutine consensus_tree(trees, p, consensus, support, info, rooted)
      !! Builds the strict or majority-rule consensus topology using ape-compatible clade-frequency selection.
      type(phylo_tree), intent(in) :: trees(:) !! Input trees with identical tip numbering and at least one tree.
      real(dp), intent(in) :: p !! Minimum clade proportion in [0.5,1]; p=0.5 uses ape's strict-majority epsilon rule.
      type(phylo_tree), intent(out) :: consensus !! Consensus topology without branch lengths.
      real(dp), allocatable, intent(out) :: support(:) !! Selected clade frequencies in consensus internal-node order.
      integer, intent(out) :: info !! Zero on success; 1 bad trees, 2 tip mismatch, 3 bad p, or 4 incompatible clades.
      logical, intent(in), optional :: rooted !! True uses rooted clades; false uses ape unrooted SHORTwise consensus.
      type(split_collection) :: collection
      logical, allocatable :: selected(:, :)
      logical, allocatable :: swap_mask(:)
      integer, allocatable :: selected_frequency(:)
      integer, allocatable :: edge(:, :)
      integer :: best_size
      integer :: i
      integer :: j
      integer :: k
      integer :: m
      integer :: n_tip
      integer :: ntree
      integer :: parent_index
      integer :: pos
      integer :: selected_count
      integer :: swap_frequency
      logical :: choose
      logical :: rooted_mode

      info = 0
      rooted_mode = .false.
      if (present(rooted)) rooted_mode = rooted
      ntree = size(trees)
      if (ntree == 0) then
         info = 1
         return
      end if
      if (p < 0.5_dp .or. p > 1.0_dp) then
         info = 3
         return
      end if
      if (.not. trees(1)%valid()) then
         info = 1
         return
      end if
      n_tip = trees(1)%n_tip
      do i = 1, ntree
         if (.not. trees(i)%valid()) then
            info = 1
            return
         end if
         if (trees(i)%n_tip /= n_tip) then
            info = 2
            return
         end if
      end do
      if (rooted_mode) then
         call prop_part(trees, collection, info)
      else
         call shortwise_part(trees, collection, info)
         if (info == 0) then
            do i = 1, collection%n_split
               if (count(collection%split(i, :)) == 0) collection%split(i, :) = .true.
            end do
         end if
      end if
      if (info /= 0) return
      selected_count = 0
      do i = 1, collection%n_split
         if (count(collection%split(i, :)) == 1) cycle
         if (p <= 0.5_dp + epsilon(p)) then
            choose = 2 * collection%frequency(i) > ntree
         else
            choose = real(collection%frequency(i), dp) >= p * real(ntree, dp)
         end if
         if (choose) selected_count = selected_count + 1
      end do
      if (selected_count == 0) then
         info = 4
         return
      end if
      allocate(selected(selected_count, n_tip), selected_frequency(selected_count))
      m = 0
      do i = 1, collection%n_split
         if (count(collection%split(i, :)) == 1) cycle
         if (p <= 0.5_dp + epsilon(p)) then
            choose = 2 * collection%frequency(i) > ntree
         else
            choose = real(collection%frequency(i), dp) >= p * real(ntree, dp)
         end if
         if (.not. choose) cycle
         m = m + 1
         selected(m, :) = collection%split(i, :)
         selected_frequency(m) = collection%frequency(i)
      end do
      allocate(swap_mask(n_tip))
      do i = 1, m - 1
         do j = i + 1, m
            if (count(selected(j, :)) > count(selected(i, :))) then
               swap_mask = selected(i, :)
               selected(i, :) = selected(j, :)
               selected(j, :) = swap_mask
               swap_frequency = selected_frequency(i)
               selected_frequency(i) = selected_frequency(j)
               selected_frequency(j) = swap_frequency
            end if
         end do
      end do
      if (.not. all(selected(1, :))) then
         info = 4
         return
      end if
      do i = 1, m - 1
         do j = i + 1, m
            if (.not. any(selected(i, :) .and. selected(j, :))) cycle
            if (all((.not. selected(j, :)) .or. selected(i, :))) cycle
            if (all((.not. selected(i, :)) .or. selected(j, :))) cycle
            info = 4
            return
         end do
      end do
      allocate(edge(n_tip + m - 1, 2))
      pos = 0
      do i = 2, m
         parent_index = 0
         best_size = n_tip + 1
         do j = 1, m
            if (count(selected(j, :)) <= count(selected(i, :))) cycle
            if (.not. all((.not. selected(i, :)) .or. selected(j, :))) cycle
            if (count(selected(j, :)) < best_size) then
               best_size = count(selected(j, :))
               parent_index = j
            end if
         end do
         if (parent_index == 0) then
            info = 4
            return
         end if
         pos = pos + 1
         edge(pos, :) = [n_tip + parent_index, n_tip + i]
      end do
      do k = 1, n_tip
         parent_index = 0
         best_size = n_tip + 1
         do i = 1, m
            if (.not. selected(i, k)) cycle
            if (count(selected(i, :)) < best_size) then
               best_size = count(selected(i, :))
               parent_index = i
            end if
         end do
         if (parent_index == 0) then
            info = 4
            return
         end if
         pos = pos + 1
         edge(pos, :) = [n_tip + parent_index, k]
      end do
      consensus = make_phylo_tree(n_tip, edge)
      if (.not. consensus%valid()) then
         info = 4
         return
      end if
      allocate(support(m))
      support = real(selected_frequency, dp) / real(ntree, dp)
   end subroutine consensus_tree

   pure subroutine tree_bipartitions(tree, splits, info)
      !! Extracts unique nontrivial edge bipartitions using ape's canonical orientation containing tip 1.
      type(phylo_tree), intent(in) :: tree !! Tree interpreted as an undirected topology while retaining its directed storage.
      logical, allocatable, intent(out) :: splits(:, :) !! Logical matrix `(nsplit,n_tip)` of canonical nontrivial bipartitions.
      integer, intent(out) :: info !! Status code: zero on success, 1 for an invalid tree.
      logical, allocatable :: descendants(:)
      logical, allocatable :: storage(:, :)
      integer, allocatable :: frequency(:)
      integer :: e
      integer :: n_side
      integer :: used

      info = 0
      if (.not. tree%valid()) then
         info = 1
         return
      end if
      allocate(storage(max(0, tree%nedge() - tree%n_tip), tree%n_tip))
      allocate(frequency(size(storage, 1)))
      storage = .false.
      frequency = 0
      used = 0
      do e = 1, tree%nedge()
         call clade_tips(tree, tree%edge(e, 2), descendants)
         n_side = count(descendants)
         if (n_side <= 1 .or. n_side >= tree%n_tip - 1) cycle
         if (.not. descendants(1)) descendants = .not. descendants
         if (.not. split_present(descendants, storage, used)) then
            used = used + 1
            storage(used, :) = descendants
         end if
      end do
      allocate(splits(used, tree%n_tip))
      if (used > 0) splits = storage(1:used, :)
   end subroutine tree_bipartitions

   pure subroutine shortwise_part(trees, collection, info)
      !! Collects rooted internal clades after ape SHORTwise canonicalization, counting each split at most once per tree.
      type(phylo_tree), intent(in) :: trees(:) !! Input trees sharing the same tip numbering.
      type(split_collection), intent(out) :: collection !! Unique SHORTwise clades and their per-tree frequencies.
      integer, intent(out) :: info !! Status code: zero on success, 1 for invalid trees, or 2 for inconsistent tip counts.
      logical, allocatable :: mask(:)
      logical, allocatable :: one_tree(:, :)
      logical, allocatable :: storage(:, :)
      integer, allocatable :: frequency(:)
      integer :: i
      integer :: j
      integer :: k
      integer :: max_split
      integer :: n_tip
      integer :: node
      integer :: used
      logical :: duplicate

      collection%n_tip = 0
      collection%n_split = 0
      info = 0
      if (size(trees) == 0) then
         allocate(collection%split(0, 0), collection%frequency(0))
         return
      end if
      if (.not. trees(1)%valid()) then
         info = 1
         return
      end if
      n_tip = trees(1)%n_tip
      max_split = 0
      do i = 1, size(trees)
         if (.not. trees(i)%valid()) then
            info = 1
            return
         end if
         if (trees(i)%n_tip /= n_tip) then
            info = 2
            return
         end if
         max_split = max_split + trees(i)%n_node
      end do
      allocate(storage(max_split, n_tip), frequency(max_split))
      storage = .false.
      frequency = 0
      used = 0
      do i = 1, size(trees)
         allocate(one_tree(trees(i)%n_node, n_tip))
         one_tree = .false.
         do j = 1, trees(i)%n_node
            node = n_tip + j
            call clade_tips(trees(i), node, mask)
            call shortwise_mask(mask)
            one_tree(j, :) = mask
         end do
         do j = 1, size(one_tree, 1)
            duplicate = .false.
            do k = 1, j - 1
               if (all(one_tree(j, :) .eqv. one_tree(k, :))) then
                  duplicate = .true.
                  exit
               end if
            end do
            if (.not. duplicate) call store_split(one_tree(j, :), storage, frequency, used)
         end do
         deallocate(one_tree)
      end do
      call finish_collection(n_tip, storage, frequency, used, collection)
   end subroutine shortwise_part

   pure subroutine shortwise_mask(mask)
      !! Canonicalizes a clade to ape's SHORTwise split representative.
      logical, intent(inout) :: mask(:) !! Tip-membership mask replaced by its complement when ape SHORTwise selects that side.
      integer :: n_member
      integer :: n_tip

      n_tip = size(mask)
      n_member = count(mask)
      if (2 * n_member > n_tip) then
         mask = .not. mask
      else if (2 * n_member == n_tip .and. n_tip > 0) then
         if (.not. mask(1)) mask = .not. mask
      end if
   end subroutine shortwise_mask

   pure subroutine store_split(mask, storage, frequency, used)
      !! Adds a split to mutable collection storage or increments its existing frequency.
      logical, intent(in) :: mask(:) !! Tip-membership mask to insert or count.
      logical, intent(inout) :: storage(:, :) !! Preallocated split storage with one split per row.
      integer, intent(inout) :: frequency(:) !! Mutable occurrence count for each occupied storage row.
      integer, intent(inout) :: used !! Number of currently occupied rows, incremented when a new split is stored.
      integer :: i

      do i = 1, used
         if (all(mask .eqv. storage(i, :))) then
            frequency(i) = frequency(i) + 1
            return
         end if
      end do
      if (used >= size(storage, 1)) return
      used = used + 1
      storage(used, :) = mask
      frequency(used) = 1
   end subroutine store_split

   pure logical function split_present(mask, storage, used) result(found)
      !! Tests whether a split mask is already present in occupied collection storage.
      logical, intent(in) :: mask(:) !! Candidate split-membership mask.
      logical, intent(in) :: storage(:, :) !! Split storage searched rowwise.
      integer, intent(in) :: used !! Number of occupied rows in `storage`.
      integer :: i

      found = .false.
      do i = 1, used
         if (all(mask .eqv. storage(i, :))) then
            found = .true.
            return
         end if
      end do
   end function split_present

   pure subroutine finish_collection(n_tip, storage, frequency, used, collection)
      !! Copies the occupied prefix of mutable storage into a compact split collection.
      integer, intent(in) :: n_tip !! Number of tips represented by every split row.
      logical, intent(in) :: storage(:, :) !! Preallocated working split storage.
      integer, intent(in) :: frequency(:) !! Working split frequencies parallel to `storage`.
      integer, intent(in) :: used !! Number of occupied rows copied to the result.
      type(split_collection), intent(out) :: collection !! Compact collection receiving unique split rows and frequencies.

      collection%n_tip = n_tip
      collection%n_split = used
      allocate(collection%split(used, n_tip), collection%frequency(used))
      if (used > 0) then
         collection%split = storage(1:used, :)
         collection%frequency = frequency(1:used)
      end if
   end subroutine finish_collection

end module ape_splits
