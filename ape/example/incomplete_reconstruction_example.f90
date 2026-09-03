program incomplete_reconstruction_example
   use ape, only : dp, phylo_tree, njs, dist_nodes
   implicit none

   real(dp) :: distance(5, 5)
   real(dp), allocatable :: fitted_distance(:, :)
   type(phylo_tree) :: tree
   integer :: info

   distance = reshape([ &
      0.0_dp, 5.0_dp, -1.0_dp, 9.0_dp, 8.0_dp, &
      5.0_dp, 0.0_dp, 10.0_dp, 10.0_dp, 9.0_dp, &
      -1.0_dp, 10.0_dp, 0.0_dp, 8.0_dp, 7.0_dp, &
      9.0_dp, 10.0_dp, 8.0_dp, 0.0_dp, 3.0_dp, &
      8.0_dp, 9.0_dp, 7.0_dp, 3.0_dp, 0.0_dp], [5, 5])

   call njs(distance, tree, info)
   if (info /= 0) error stop 'NJ* reconstruction failed'

   call dist_nodes(tree, fitted_distance, info)
   if (info /= 0) error stop 'patristic-distance calculation failed'

   print '(a,f6.2)', 'NJ* inferred the missing distance d(1,3) = ', fitted_distance(1, 3)
end program incomplete_reconstruction_example
