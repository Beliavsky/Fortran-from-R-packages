program reconstruct_example
   use ape
   implicit none
   type(phylo_tree) :: tree
   type(reconstruct_result) :: fit
   integer :: edge(6, 2)
   real(dp) :: x(4)
   integer :: info

   edge = reshape([5, 5, 6, 6, 7, 7, 6, 7, 1, 2, 3, 4], [6, 2])
   tree = make_phylo_tree(4, edge, [1.0_dp, 0.8_dp, 0.7_dp, 1.2_dp, 1.0_dp, 1.4_dp])
   x = [1.0_dp, 3.0_dp, 5.0_dp, 7.0_dp]

   call reconstruct_fit(x, tree, 'GLS_OU', fit, info, alpha=0.4_dp)
   if (info /= 0) error stop 'reconstruct example failed'
   print '(a,*(f10.5,1x))', 'ancestral estimates: ', fit%ancestral
   print '(a,f10.5)', 'OU theta: ', fit%theta
   print '(a,f12.6)', 'log likelihood: ', fit%log_likelihood
end program reconstruct_example
