program discrete_ace_example
   use ape
   implicit none
   type(phylo_tree) :: tree
   type(ace_discrete_result) :: fit
   integer :: edge(6, 2)
   integer :: info

   edge = reshape([5, 5, 6, 6, 7, 7, 6, 7, 1, 2, 3, 4], [6, 2])
   tree = make_phylo_tree(4, edge, [1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp])
   call ace_discrete_fit(tree, [1, 1, 2, 2], 2, 'ER', fit, info)
   if (info /= 0) error stop 'discrete ACE fit failed'
   print '(a,f10.6)', 'ER transition rate: ', fit%rates(1)
   print '(a,f10.6)', 'log likelihood: ', fit%log_likelihood
   print '(a,2f10.6)', 'root state likelihoods: ', fit%ancestral_likelihood(1, :)
end program discrete_ace_example
