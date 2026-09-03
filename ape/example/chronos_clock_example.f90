program chronos_clock_example
   use ape
   implicit none
   type(phylo_tree) :: tree
   type(chronos_clock_result) :: fit
   integer :: edge(6, 2)
   integer :: info

   edge = reshape([5, 5, 6, 6, 7, 7, 6, 7, 1, 2, 3, 4], [6, 2])
   tree = make_phylo_tree(4, edge, [0.25_dp, 0.25_dp, 0.25_dp, 0.25_dp, 0.25_dp, 0.25_dp])
   call chronos_clock_fit(tree, fit, info)
   if (info /= 0) error stop 'chronos clock fit failed'
   print '(a,f10.6)', 'common clock rate: ', fit%rate
   print '(a,3f10.6)', 'internal ages: ', fit%node_age
   print '(a,f10.6)', 'log likelihood: ', fit%log_likelihood
   print '(a,f10.6)', 'PHIIC: ', fit%phiic
end program chronos_clock_example
