! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
!
! Tree-distance/depth algorithms are derived from computational routines in
! ape, including src/dist_nodes.c, src/plot_phylo.c, R/branching.times.R,
! R/nodepath.R, R/balance.R, and R/cherry.R.
! Independent contrasts follow src/pic.c (Copyright 2006-2017 Emmanuel Paradis).
module ape_tree_algorithms
   use r_kinds, only : dp
   use ape_types, only : phylo_tree, parent_vector, child_counts, edge_index_to_child
   implicit none
   private

   public :: node_depth_edgelength
   public :: node_depth_count
   public :: dist_nodes
   public :: mrca
   public :: node_path
   public :: descendant_tip_counts
   public :: balance_counts
   public :: cherry_count
   public :: branching_times
   public :: pic
   public :: ace_pic
   public :: chrono_mpl
   public :: compute_brtime

contains

   pure subroutine node_depth_edgelength(tree, depth, info)
      !! Computes root-to-node path lengths for all tips and internal nodes.
      type(phylo_tree), intent(in) :: tree !! Rooted phylogenetic tree with one branch length per edge.
      real(dp), allocatable, intent(out) :: depth(:) !! Root-to-node distances indexed by ape node number.
      integer, intent(out) :: info !! Status code: zero on success, nonzero for invalid tree or missing lengths.
      integer, allocatable :: edge_for_child(:)
      integer, allocatable :: parent(:)
      integer :: child
      integer :: e
      integer :: node
      integer :: n_total
      integer :: steps
      real(dp) :: value

      info = 0
      n_total = tree%total_nodes()
      allocate(depth(n_total))
      depth = 0.0_dp
      if (.not. tree%valid()) then
         info = 1
         return
      end if
      if (.not. tree%has_lengths()) then
         info = 2
         return
      end if
      parent = parent_vector(tree)
      edge_for_child = edge_index_to_child(tree)
      do node = 1, n_total
         child = node
         value = 0.0_dp
         steps = 0
         do while (parent(child) /= 0)
            e = edge_for_child(child)
            if (e == 0) then
               info = 3
               return
            end if
            value = value + tree%edge_length(e)
            child = parent(child)
            steps = steps + 1
            if (steps > n_total) then
               info = 4
               return
            end if
         end do
         depth(node) = value
      end do
   end subroutine node_depth_edgelength

   pure subroutine descendant_tip_counts(tree, counts, info)
      !! Counts how many terminal taxa descend from every node.
      type(phylo_tree), intent(in) :: tree !! Rooted tree whose descendant-tip counts are required.
      integer, allocatable, intent(out) :: counts(:) !! Number of descendant tips for each ape node number.
      integer, intent(out) :: info !! Status code: zero on success, nonzero if the tree encoding is invalid.
      integer, allocatable :: parent(:)
      integer :: node
      integer :: current
      integer :: n_total
      integer :: steps

      info = 0
      n_total = tree%total_nodes()
      allocate(counts(n_total))
      counts = 0
      if (.not. tree%valid()) then
         info = 1
         return
      end if
      parent = parent_vector(tree)
      do node = 1, tree%n_tip
         current = node
         steps = 0
         do
            counts(current) = counts(current) + 1
            if (parent(current) == 0) exit
            current = parent(current)
            steps = steps + 1
            if (steps > n_total) then
               info = 2
               return
            end if
         end do
      end do
   end subroutine descendant_tip_counts

   pure subroutine node_depth_count(tree, method, depth, info)
      !! Reproduces ape node.depth topology-based depths for methods 1 and 2.
      type(phylo_tree), intent(in) :: tree !! Rooted tree whose topology-based depths are required.
      integer, intent(in) :: method !! ape method: 1 counts descendant tips; 2 is maximum descendant edge depth plus one.
      integer, allocatable, intent(out) :: depth(:) !! Integer depth values indexed by ape node number.
      integer, intent(out) :: info !! Status code: zero on success, nonzero for an invalid tree or method.
      integer, allocatable :: child_count(:)
      integer, allocatable :: parent(:)
      logical, allocatable :: done(:)
      integer :: e
      integer :: i
      integer :: n_total
      integer :: processed
      logical :: progressed

      info = 0
      if (method == 1) then
         call descendant_tip_counts(tree, depth, info)
         return
      end if
      n_total = tree%total_nodes()
      allocate(depth(n_total), done(n_total))
      depth = 0
      done = .false.
      if (.not. tree%valid()) then
         info = 1
         return
      end if
      if (method /= 2) then
         info = 2
         return
      end if
      child_count = child_counts(tree)
      parent = parent_vector(tree)
      do i = 1, tree%n_tip
         depth(i) = 1
         done(i) = .true.
      end do
      processed = tree%n_tip
      do while (processed < n_total)
         progressed = .false.
         do i = tree%n_tip + 1, n_total
            if (done(i)) cycle
            if (all_children_done(tree, i, done)) then
               do e = 1, tree%nedge()
                  if (tree%edge(e, 1) == i) depth(i) = max(depth(i), depth(tree%edge(e, 2)) + 1)
               end do
               if (child_count(i) == 0) depth(i) = 0
               done(i) = .true.
               processed = processed + 1
               progressed = .true.
            end if
         end do
         if (.not. progressed) then
            info = 3
            return
         end if
      end do
   end subroutine node_depth_count

   pure elemental integer function mrca(tree, node_a, node_b) result(ancestor)
      !! Finds the most recent common ancestor of two nodes in a rooted tree.
      type(phylo_tree), intent(in) :: tree !! Rooted tree containing both requested nodes.
      integer, intent(in) :: node_a !! First ape node number, from 1 through the total number of nodes.
      integer, intent(in) :: node_b !! Second ape node number, from 1 through the total number of nodes.
      integer, allocatable :: parent(:)
      logical, allocatable :: lineage(:)
      integer :: current
      integer :: n_total
      integer :: steps

      ancestor = 0
      n_total = tree%total_nodes()
      if (.not. tree%valid()) return
      if (node_a < 1 .or. node_a > n_total) return
      if (node_b < 1 .or. node_b > n_total) return
      parent = parent_vector(tree)
      allocate(lineage(n_total))
      lineage = .false.
      current = node_a
      steps = 0
      do
         lineage(current) = .true.
         if (parent(current) == 0) exit
         current = parent(current)
         steps = steps + 1
         if (steps > n_total) return
      end do
      current = node_b
      steps = 0
      do
         if (lineage(current)) then
            ancestor = current
            return
         end if
         if (parent(current) == 0) exit
         current = parent(current)
         steps = steps + 1
         if (steps > n_total) return
      end do
   end function mrca

   pure subroutine node_path(tree, from_node, to_node, path, info)
      !! Returns the inclusive simple path between two nodes.
      type(phylo_tree), intent(in) :: tree !! Rooted tree containing the requested endpoints.
      integer, intent(in) :: from_node !! Starting ape node number.
      integer, intent(in) :: to_node !! Ending ape node number.
      integer, allocatable, intent(out) :: path(:) !! Inclusive node sequence from `from_node` to `to_node`.
      integer, intent(out) :: info !! Status code: zero on success, nonzero for invalid nodes/tree.
      integer, allocatable :: parent(:)
      integer, allocatable :: left(:)
      integer, allocatable :: right(:)
      integer :: ancestor
      integer :: current
      integer :: i
      integer :: nl
      integer :: nr
      integer :: n_total

      info = 0
      n_total = tree%total_nodes()
      if (.not. tree%valid()) then
         allocate(path(0))
         info = 1
         return
      end if
      if (from_node < 1 .or. from_node > n_total .or. to_node < 1 .or. to_node > n_total) then
         allocate(path(0))
         info = 2
         return
      end if
      ancestor = mrca(tree, from_node, to_node)
      if (ancestor == 0) then
         allocate(path(0))
         info = 3
         return
      end if
      parent = parent_vector(tree)
      allocate(left(n_total), right(n_total))
      nl = 0
      current = from_node
      do
         nl = nl + 1
         left(nl) = current
         if (current == ancestor) exit
         current = parent(current)
      end do
      nr = 0
      current = to_node
      do while (current /= ancestor)
         nr = nr + 1
         right(nr) = current
         current = parent(current)
      end do
      allocate(path(nl + nr))
      path(1:nl) = left(1:nl)
      do i = 1, nr
         path(nl + i) = right(nr - i + 1)
      end do
   end subroutine node_path

   pure subroutine dist_nodes(tree, distance, info)
      !! Computes all pairwise patristic distances among tips and internal nodes.
      type(phylo_tree), intent(in) :: tree !! Rooted tree with branch lengths.
      real(dp), allocatable, intent(out) :: distance(:, :) !! Symmetric pairwise node-distance matrix.
      integer, intent(out) :: info !! Status code: zero on success, nonzero for invalid tree or branch lengths.
      real(dp), allocatable :: depth(:)
      integer :: ancestor
      integer :: i
      integer :: j
      integer :: n_total

      call node_depth_edgelength(tree, depth, info)
      n_total = tree%total_nodes()
      allocate(distance(n_total, n_total))
      distance = 0.0_dp
      if (info /= 0) return
      do i = 1, n_total - 1
         do j = i + 1, n_total
            ancestor = mrca(tree, i, j)
            if (ancestor == 0) then
               info = 4
               return
            end if
            distance(i, j) = depth(i) + depth(j) - 2.0_dp * depth(ancestor)
            distance(j, i) = distance(i, j)
         end do
      end do
   end subroutine dist_nodes

   pure subroutine balance_counts(tree, balances, info)
      !! Returns descendant-tip counts for the two children of each internal node.
      type(phylo_tree), intent(in) :: tree !! Rooted binary tree whose node balances are requested.
      integer, allocatable, intent(out) :: balances(:, :) !! Shape `(n_node, 2)`; rows follow node numbers `n_tip+1:`.
      integer, intent(out) :: info !! Status code: zero on success; nonzero if the tree is invalid or nonbinary.
      integer, allocatable :: counts(:)
      integer, allocatable :: nchild(:)
      integer :: e
      integer :: node
      integer :: slot

      call descendant_tip_counts(tree, counts, info)
      allocate(balances(tree%n_node, 2))
      balances = 0
      if (info /= 0) return
      nchild = child_counts(tree)
      do node = tree%n_tip + 1, tree%total_nodes()
         if (nchild(node) /= 2) then
            info = 2
            return
         end if
         slot = 0
         do e = 1, tree%nedge()
            if (tree%edge(e, 1) /= node) cycle
            slot = slot + 1
            balances(node - tree%n_tip, slot) = counts(tree%edge(e, 2))
         end do
      end do
   end subroutine balance_counts

   pure elemental integer function cherry_count(tree) result(count)
      !! Counts binary internal nodes whose two direct descendants are both tips.
      type(phylo_tree), intent(in) :: tree !! Rooted tree in ape node numbering.
      integer, allocatable :: tip_children(:)
      integer, allocatable :: total_children(:)
      integer :: e
      integer :: parent

      count = 0
      if (.not. tree%valid()) return
      allocate(tip_children(tree%total_nodes()))
      tip_children = 0
      total_children = child_counts(tree)
      do e = 1, tree%nedge()
         parent = tree%edge(e, 1)
         if (tree%edge(e, 2) <= tree%n_tip) tip_children(parent) = tip_children(parent) + 1
      end do
      do parent = tree%n_tip + 1, tree%total_nodes()
         if (total_children(parent) == 2 .and. tip_children(parent) == 2) count = count + 1
      end do
   end function cherry_count

   pure subroutine branching_times(tree, times, info)
      !! Computes ape-style ages of internal nodes measured backward from the most distant tip.
      type(phylo_tree), intent(in) :: tree !! Rooted tree with branch lengths, typically ultrametric for biological ages.
      real(dp), allocatable, intent(out) :: times(:) !! Ages for internal nodes `n_tip+1:` in node-number order.
      integer, intent(out) :: info !! Status code: zero on success, nonzero if depths cannot be computed.
      real(dp), allocatable :: depth(:)
      real(dp) :: tree_depth
      integer :: i

      call node_depth_edgelength(tree, depth, info)
      allocate(times(tree%n_node))
      times = 0.0_dp
      if (info /= 0) return
      tree_depth = maxval(depth(1:tree%n_tip))
      do i = 1, tree%n_node
         times(i) = tree_depth - depth(tree%n_tip + i)
      end do
   end subroutine branching_times

   pure subroutine pic(tree, phenotype, contrast, variance, scaled, rescaled_tree, info, ancestral)
      !! Computes Felsenstein phylogenetically independent contrasts on a rooted binary tree.
      type(phylo_tree), intent(in) :: tree !! Rooted binary tree with positive branch lengths.
      real(dp), intent(in) :: phenotype(:) !! Tip trait values, one value for each tip in node-number order.
      real(dp), allocatable, intent(out) :: contrast(:) !! One contrast per internal node, indexed by internal node number.
      real(dp), allocatable, intent(out) :: variance(:) !! Sum of the two daughter branch lengths for each contrast.
      logical, intent(in), optional :: scaled !! If true, divide each contrast by the square root of its variance.
      type(phylo_tree), intent(out), optional :: rescaled_tree !! Tree with ancestral-edge lengths updated as in ape `pic`.
      integer, intent(out) :: info !! Status code: zero on success, nonzero for invalid input or nonbinary topology.
      real(dp), allocatable, intent(out), optional :: ancestral(:) !! Estimated trait values at internal nodes in ape node order.
      integer, allocatable :: edge_for_child(:)
      integer, allocatable :: nchild(:)
      integer, allocatable :: child_one(:)
      integer, allocatable :: child_two(:)
      logical, allocatable :: done(:)
      real(dp), allocatable :: value(:)
      real(dp), allocatable :: edge_length(:)
      integer :: anc
      integer :: child
      integer :: e
      integer :: ic
      integer :: n_total
      integer :: processed
      real(dp) :: l1
      real(dp) :: l2
      real(dp) :: sumbl
      logical :: do_scale
      logical :: progressed

      info = 0
      allocate(contrast(tree%n_node), variance(tree%n_node))
      contrast = 0.0_dp
      variance = 0.0_dp
      if (.not. tree%valid() .or. .not. tree%has_lengths()) then
         info = 1
         return
      end if
      if (size(phenotype) /= tree%n_tip) then
         info = 2
         return
      end if
      n_total = tree%total_nodes()
      nchild = child_counts(tree)
      do anc = tree%n_tip + 1, n_total
         if (nchild(anc) /= 2) then
            info = 3
            return
         end if
      end do
      allocate(child_one(n_total), child_two(n_total), done(n_total), value(n_total))
      child_one = 0
      child_two = 0
      done = .false.
      value = 0.0_dp
      value(1:tree%n_tip) = phenotype
      done(1:tree%n_tip) = .true.
      edge_length = tree%edge_length
      edge_for_child = edge_index_to_child(tree)
      do e = 1, tree%nedge()
         anc = tree%edge(e, 1)
         child = tree%edge(e, 2)
         if (child_one(anc) == 0) then
            child_one(anc) = child
         else
            child_two(anc) = child
         end if
      end do
      do_scale = .true.
      if (present(scaled)) do_scale = scaled
      processed = 0
      do while (processed < tree%n_node)
         progressed = .false.
         do anc = tree%n_tip + 1, n_total
            if (done(anc)) cycle
            if (.not. done(child_one(anc)) .or. .not. done(child_two(anc))) cycle
            l1 = edge_length(edge_for_child(child_one(anc)))
            l2 = edge_length(edge_for_child(child_two(anc)))
            sumbl = l1 + l2
            if (sumbl <= 0.0_dp) then
               info = 4
               return
            end if
            ic = anc - tree%n_tip
            contrast(ic) = value(child_one(anc)) - value(child_two(anc))
            variance(ic) = sumbl
            if (do_scale) contrast(ic) = contrast(ic) / sqrt(sumbl)
            value(anc) = (value(child_one(anc)) * l2 + value(child_two(anc)) * l1) / sumbl
            if (edge_for_child(anc) /= 0) then
               e = edge_for_child(anc)
               edge_length(e) = edge_length(e) + l1 * l2 / sumbl
            end if
            done(anc) = .true.
            processed = processed + 1
            progressed = .true.
         end do
         if (.not. progressed) then
            info = 5
            return
         end if
      end do
      if (present(rescaled_tree)) then
         rescaled_tree = tree
         rescaled_tree%edge_length = edge_length
      end if
      if (present(ancestral)) then
         allocate(ancestral(tree%n_node))
         ancestral = value(tree%n_tip + 1:n_total)
      end if
   end subroutine pic

   pure subroutine ace_pic(tree, phenotype, estimates, variance, scaled, info)
      !! Estimates continuous ancestral states with ape's Brownian-motion PIC method.
      type(phylo_tree), intent(in) :: tree !! Rooted binary tree with positive branch lengths.
      real(dp), intent(in) :: phenotype(:) !! Tip trait values in node-number order, with one value per tip.
      real(dp), allocatable, intent(out) :: estimates(:) !! Internal-node ancestral estimates in ape node-number order.
      real(dp), allocatable, intent(out) :: variance(:) !! Contrast variances used by ape for PIC confidence intervals.
      logical, intent(in), optional :: scaled !! If true, compute scaled contrasts internally; ancestral estimates are unchanged.
      integer, intent(out) :: info !! Status code: zero on success, otherwise the corresponding `pic` validation failure.
      real(dp), allocatable :: contrast(:)

      if (present(scaled)) then
         call pic(tree, phenotype, contrast, variance, scaled=scaled, info=info, ancestral=estimates)
      else
         call pic(tree, phenotype, contrast, variance, info=info, ancestral=estimates)
      end if
   end subroutine ace_pic

   pure subroutine chrono_mpl(tree, dated_tree, info, standard_error, p_value)
      !! Dates a binary tree by ape's mean-path-length method `chronoMPL`.
      type(phylo_tree), intent(in) :: tree !! Rooted binary tree with branch lengths interpreted as substitution path lengths.
      type(phylo_tree), intent(out) :: dated_tree !! Copy whose edge lengths are parent-minus-child mean path ages.
      integer, intent(out) :: info !! Status: zero on success, 1 invalid/missing lengths, 2 nonbinary, or 3 cyclic/unresolved order.
      real(dp), allocatable, intent(out), optional :: standard_error(:) !! Internal-node age SEs in ape node-number order.
      real(dp), allocatable, intent(out), optional :: p_value(:) !! Internal-node two-sided molecular-clock tests from `chronoMPL`.
      integer, allocatable :: child_one(:)
      integer, allocatable :: child_two(:)
      integer, allocatable :: edge_for_child(:)
      integer, allocatable :: ndesc(:)
      logical, allocatable :: done(:)
      real(dp), allocatable :: path_sum(:)
      real(dp), allocatable :: path_ss(:)
      real(dp), allocatable :: node_age(:)
      real(dp), allocatable :: local_p(:)
      integer :: anc
      integer :: e
      integer :: n_total
      integer :: processed
      integer :: d1
      integer :: d2
      real(dp) :: branch_a
      real(dp) :: branch_b
      real(dp) :: a
      real(dp) :: b
      real(dp) :: z
      real(dp) :: z_variance
      logical :: progressed

      dated_tree = tree
      info = 0
      if (.not. tree%valid() .or. .not. tree%has_lengths()) then
         info = 1
         return
      end if
      n_total = tree%total_nodes()
      ndesc = child_counts(tree)
      do anc = tree%n_tip + 1, n_total
         if (ndesc(anc) /= 2) then
            info = 2
            return
         end if
      end do
      call descendant_tip_counts(tree, ndesc, info)
      if (info /= 0) then
         info = 3
         return
      end if
      allocate(child_one(n_total), child_two(n_total), done(n_total))
      allocate(path_sum(n_total), path_ss(n_total), node_age(n_total), local_p(tree%n_node))
      child_one = 0
      child_two = 0
      done = .false.
      done(1:tree%n_tip) = .true.
      path_sum = 0.0_dp
      path_ss = 0.0_dp
      node_age = 0.0_dp
      local_p = 0.0_dp
      edge_for_child = edge_index_to_child(tree)
      do e = 1, tree%nedge()
         anc = tree%edge(e, 1)
         if (child_one(anc) == 0) then
            child_one(anc) = tree%edge(e, 2)
         else
            child_two(anc) = tree%edge(e, 2)
         end if
      end do

      processed = 0
      do while (processed < tree%n_node)
         progressed = .false.
         do anc = tree%n_tip + 1, n_total
            if (done(anc)) cycle
            d1 = child_one(anc)
            d2 = child_two(anc)
            if (.not. done(d1) .or. .not. done(d2)) cycle
            branch_a = tree%edge_length(edge_for_child(d1))
            branch_b = tree%edge_length(edge_for_child(d2))
            a = path_sum(d1) + real(ndesc(d1), dp) * branch_a
            b = path_sum(d2) + real(ndesc(d2), dp) * branch_b
            path_sum(anc) = a + b
            path_ss(anc) = path_ss(d1) + real(ndesc(d1) * ndesc(d1), dp) * branch_a + &
               path_ss(d2) + real(ndesc(d2) * ndesc(d2), dp) * branch_b
            z = abs(a / real(ndesc(d1), dp) - b / real(ndesc(d2), dp))
            z_variance = (path_ss(d1) + real(ndesc(d1) * ndesc(d1), dp) * branch_a) / &
               real(ndesc(d1) * ndesc(d1), dp)
            z_variance = z_variance + &
               (path_ss(d2) + real(ndesc(d2) * ndesc(d2), dp) * branch_b) / &
               real(ndesc(d2) * ndesc(d2), dp)
            if (z_variance > 0.0_dp) then
               local_p(anc - tree%n_tip) = erfc((z / sqrt(z_variance)) / sqrt(2.0_dp))
            else if (abs(z) <= tiny(1.0_dp)) then
               local_p(anc - tree%n_tip) = 1.0_dp
            else
               local_p(anc - tree%n_tip) = 0.0_dp
            end if
            node_age(anc) = path_sum(anc) / real(ndesc(anc), dp)
            done(anc) = .true.
            processed = processed + 1
            progressed = .true.
         end do
         if (.not. progressed) then
            info = 3
            return
         end if
      end do

      do e = 1, tree%nedge()
         dated_tree%edge_length(e) = node_age(tree%edge(e, 1)) - node_age(tree%edge(e, 2))
      end do
      if (present(standard_error)) then
         allocate(standard_error(tree%n_node))
         do anc = tree%n_tip + 1, n_total
            standard_error(anc - tree%n_tip) = sqrt(path_ss(anc)) / real(ndesc(anc), dp)
         end do
      end if
      if (present(p_value)) then
         allocate(p_value(tree%n_node))
         p_value = local_p
      end if
   end subroutine chrono_mpl

   pure subroutine compute_brtime(tree, branching_age, dated_tree, info)
      !! Applies numeric node ages as ape `compute.brtime(phy, method=<numeric>, force.positive=FALSE)`.
      type(phylo_tree), intent(in) :: tree !! Tree topology whose edge lengths are replaced; existing lengths are not required.
      real(dp), intent(in) :: branching_age(:) !! Ages for internal nodes `n_tip+1:` in node-number order; tips have age zero.
      type(phylo_tree), intent(out) :: dated_tree !! Tree copy with each edge length equal to parent age minus child age.
      integer, intent(out) :: info !! Status: zero on success, 1 invalid tree, or 2 if the number of ages differs from `n_node`.
      real(dp), allocatable :: node_age(:)
      integer :: e

      dated_tree = tree
      info = 0
      if (.not. tree%valid()) then
         info = 1
         return
      end if
      if (size(branching_age) /= tree%n_node) then
         info = 2
         return
      end if
      allocate(node_age(tree%total_nodes()))
      node_age = 0.0_dp
      node_age(tree%n_tip + 1:) = branching_age
      if (.not. allocated(dated_tree%edge_length)) allocate(dated_tree%edge_length(tree%nedge()))
      do e = 1, tree%nedge()
         dated_tree%edge_length(e) = node_age(tree%edge(e, 1)) - node_age(tree%edge(e, 2))
      end do
   end subroutine compute_brtime

   pure logical function all_children_done(tree, node, done) result(ok)
      type(phylo_tree), intent(in) :: tree !! Tree containing the child edges to inspect.
      integer, intent(in) :: node !! Parent ape node number whose children are tested.
      logical, intent(in) :: done(:) !! Node-indexed mask marking already evaluated descendants.
      integer :: e

      ok = .true.
      do e = 1, tree%nedge()
         if (tree%edge(e, 1) == node) then
            if (.not. done(tree%edge(e, 2))) then
               ok = .false.
               return
            end if
         end if
      end do
   end function all_children_done

end module ape_tree_algorithms
