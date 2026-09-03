program chronopl_example
   use ape
   implicit none
   type(phylo_tree) :: tree
   type(chronopl_result) :: fit
   integer :: edge(6, 2)
   integer :: info

   edge = reshape([5, 5, 6, 6, 7, 7, 6, 7, 1, 2, 3, 4], [6, 2])
   tree = make_phylo_tree(4, edge, [0.25_dp, 0.25_dp, 0.25_dp, 0.25_dp, 0.25_dp, 0.25_dp])
   call chronopl_fit(tree, 1.0_dp, fit, info)
   if (info /= 0) error stop 'chronopl fit failed'
   print '(a,f10.6)', 'penalized log likelihood: ', fit%penalized_log_likelihood
   print '(a,3f10.6)', 'internal node ages: ', fit%node_age
   print '(a,6f10.6)', 'edge rates: ', fit%rates
end program chronopl_example
