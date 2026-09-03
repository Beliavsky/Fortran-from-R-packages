program binary_pglmm_example
   use ape
   implicit none
   type(phylo_tree) :: tree
   type(binary_pglmm_result) :: fit
   integer :: edge(14, 2)
   integer :: y(8)
   real(dp) :: design(8, 2)
   integer :: info

   edge = reshape([ &
      9, 9, 10, 10, 11, 11, 12, 12, 13, 13, 14, 14, 15, 15, &
      10, 11, 12, 13, 14, 15, 1, 2, 3, 4, 5, 6, 7, 8], [14, 2])
   tree = make_phylo_tree(8, edge, [ &
      1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, &
      1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp])
   design(:, 1) = 1.0_dp
   design(:, 2) = [-1.5_dp, -1.0_dp, -0.5_dp, 0.0_dp, 0.2_dp, 0.7_dp, 1.1_dp, 1.5_dp]
   y = [0, 1, 0, 0, 1, 1, 0, 1]

   call binary_pglmm_fit(y, design, tree, fit, info)
   if (info /= 0) error stop 'binaryPGLMM example failed'
   print '(a,f10.6)', 'phylogenetic variance s2: ', fit%s2
   print '(a,*(f10.6,1x))', 'fixed effects: ', fit%coefficients
   print '(a,f12.6)', 'conditional REML log likelihood: ', fit%conditional_reml_log_likelihood
end program binary_pglmm_example
