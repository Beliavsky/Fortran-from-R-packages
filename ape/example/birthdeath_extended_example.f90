program birthdeath_extended_example
   use ape
   implicit none
   type(phylo_tree) :: tree
   type(birthdeath_extended_result) :: fit
   integer :: edge(6, 2)
   integer :: info

   edge = reshape([5, 5, 6, 6, 7, 7, 1, 6, 2, 7, 3, 4], [6, 2])
   tree = make_phylo_tree(4, edge, [3.0_dp, 1.0_dp, 2.0_dp, 1.0_dp, 1.0_dp, 1.0_dp])
   call birthdeath_extended_fit(tree, [5, 2, 3, 4], fit, info)
   if (info /= 0) error stop 'extended birth-death fit failed'

   print '(a,f10.6)', 'd/b = ', fit%death_birth_ratio
   print '(a,f10.6)', 'b-d = ', fit%net_diversification
   print '(a,f12.6)', 'log likelihood = ', fit%log_likelihood
end program birthdeath_extended_example
