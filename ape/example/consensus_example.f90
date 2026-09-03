! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
program consensus_example
   use ape, only : dp, phylo_tree, make_phylo_tree, consensus_tree, tree_bipartitions
   implicit none

   type(phylo_tree) :: tree_a
   type(phylo_tree) :: tree_b
   type(phylo_tree) :: trees(3)
   type(phylo_tree) :: consensus
   integer :: edge_a(5, 2)
   integer :: edge_b(5, 2)
   logical, allocatable :: splits(:, :)
   real(dp), allocatable :: support(:)
   integer :: info

   edge_a = reshape([5, 5, 5, 6, 6, 1, 2, 6, 3, 4], [5, 2])
   edge_b = reshape([5, 5, 5, 6, 6, 1, 3, 6, 2, 4], [5, 2])
   tree_a = make_phylo_tree(4, edge_a)
   tree_b = make_phylo_tree(4, edge_b)
   trees = [tree_a, tree_a, tree_b]

   call consensus_tree(trees, 0.5_dp, consensus, support, info)
   if (info /= 0) error stop 'consensus_tree failed'
   call tree_bipartitions(consensus, splits, info)
   if (info /= 0) error stop 'tree_bipartitions failed'

   print '(a,i0)', 'majority consensus splits: ', size(splits, 1)
   print '(a,*(f8.4,1x))', 'node support: ', support
end program consensus_example
