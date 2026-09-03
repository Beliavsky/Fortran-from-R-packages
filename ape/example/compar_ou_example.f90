program compar_ou_example
   use ape
   implicit none
   type(phylo_tree) :: tree
   type(compar_ou_result) :: fit
   integer :: edge(6, 2)
   integer :: info

   edge = reshape([5, 5, 6, 6, 7, 7, 6, 7, 1, 2, 3, 4], [6, 2])
   tree = make_phylo_tree(4, edge, [1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp])
   call compar_ou_fit(tree, [-0.989121350348_dp, -0.367786651468_dp, 1.287925261289_dp, 0.193974419133_dp], &
      fit, info)
   if (info /= 0) error stop 'compar.ou fit failed'
   print '(a,f10.6)', 'alpha: ', fit%alpha
   print '(a,f10.6)', 'sigma2: ', fit%sigma2
   print '(a,f10.6)', 'theta: ', fit%theta(1)
   print '(a,f10.6)', 'log likelihood: ', fit%log_likelihood
end program compar_ou_example
