program ace_likelihood_example
   use ape, only : ace_continuous_ml, ace_continuous_reml, ace_continuous_result, dp, &
      make_phylo_tree, phylo_tree
   implicit none

   integer :: edge(6, 2)
   real(dp) :: edge_length(6)
   real(dp) :: phenotype(4)
   type(ace_continuous_result) :: fit
   type(phylo_tree) :: tree
   integer :: info

   edge = reshape([5, 5, 6, 6, 7, 7, 6, 7, 1, 2, 3, 4], [6, 2])
   edge_length = 1.0_dp
   phenotype = [1.0_dp, 3.0_dp, 5.0_dp, 7.0_dp]
   tree = make_phylo_tree(4, edge, edge_length)

   call ace_continuous_ml(tree, phenotype, fit, info)
   if (info /= 0) error stop 'continuous ACE ML failed'
   print '(a,f10.6)', 'ML sigma^2: ', fit%sigma2
   print '(a,*(1x,f10.6))', 'ML ancestral states:', fit%estimates

   call ace_continuous_reml(tree, phenotype, fit, info)
   if (info /= 0) error stop 'continuous ACE REML failed'
   print '(a,f10.6)', 'REML sigma^2: ', fit%sigma2
   print '(a,*(1x,f10.6))', 'REML ancestral states:', fit%estimates
end program ace_likelihood_example
