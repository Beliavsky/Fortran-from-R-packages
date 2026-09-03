program corphylo_example
   use ape
   implicit none
   type(phylo_tree) :: tree
   type(corphylo_result) :: fit
   integer :: edge(6, 2)
   real(dp) :: x(4, 2)
   integer :: info

   edge = reshape([5, 5, 6, 6, 7, 7, 6, 7, 1, 2, 3, 4], [6, 2])
   tree = make_phylo_tree(4, edge, [1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp])
   x(:, 1) = [1.0_dp, 2.0_dp, 4.0_dp, 6.0_dp]
   x(:, 2) = [2.0_dp, 1.0_dp, 5.0_dp, 4.0_dp]

   call corphylo_fit(x, tree, fit, info, constrain_d=.true.)
   if (info /= 0) error stop 'corphylo example failed'
   print '(a,f12.6)', 'log likelihood: ', fit%log_likelihood
   print '(a,*(f10.5,1x))', 'OU d: ', fit%d
   print '(a,f10.5)', 'trait correlation: ', fit%correlation(1, 2)
end program corphylo_example
