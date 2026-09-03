! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
program tree_statistics_example
   use ape, only : dp, phylo_tree, make_phylo_tree, gamma_stat, branching_times
   implicit none

   integer :: edge(6, 2)
   real(dp) :: edge_length(6)
   real(dp), allocatable :: times(:)
   real(dp) :: gamma
   type(phylo_tree) :: tree
   integer :: info

   edge = reshape([5, 5, 6, 6, 7, 7, 6, 7, 1, 2, 3, 4], [6, 2])
   edge_length = 1.0_dp
   tree = make_phylo_tree(4, edge, edge_length)

   call branching_times(tree, times, info)
   if (info /= 0) error stop 'branching_times failed'
   call gamma_stat(tree, gamma, info)
   if (info /= 0) error stop 'gamma_stat failed'

   print '(a,*(f8.3,1x))', 'branching times: ', times
   print '(a,f10.6)', 'gamma: ', gamma
end program tree_statistics_example
