program pgls_example
   use ape
   implicit none
   type(phylo_tree) :: tree
   type(pgls_result) :: fit
   real(dp) :: design(4, 2)
   integer :: edge(6, 2)
   integer :: info

   edge = reshape([5, 5, 6, 6, 7, 7, 6, 7, 1, 2, 3, 4], [6, 2])
   tree = make_phylo_tree(4, edge, [1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp])
   design(:, 1) = 1.0_dp
   design(:, 2) = [0.0_dp, 1.0_dp, 2.0_dp, 3.0_dp]
   call pgls_fit_model(tree, [1.0_dp, 2.0_dp, 4.0_dp, 8.0_dp], design, 'pagel', fit, info)
   if (info /= 0) error stop 'PGLS fit failed'
   print '(a,f10.6)', 'Pagel lambda: ', fit%correlation_parameter
   print '(a,2f10.6)', 'GLS coefficients: ', fit%coefficients
   print '(a,f10.6)', 'profile log likelihood: ', fit%log_likelihood
end program pgls_example
