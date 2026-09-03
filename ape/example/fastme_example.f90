program fastme_example
   use ape, only : dp, fastme_bal, phylo_tree
   implicit none

   real(dp) :: distance(6, 6)
   type(phylo_tree) :: tree
   integer :: info
   integer :: e

   distance = reshape([ &
      0.0_dp, 5.0_dp, 9.0_dp, 9.0_dp, 8.0_dp, 7.0_dp, &
      5.0_dp, 0.0_dp, 10.0_dp, 10.0_dp, 9.0_dp, 8.0_dp, &
      9.0_dp, 10.0_dp, 0.0_dp, 8.0_dp, 7.0_dp, 9.0_dp, &
      9.0_dp, 10.0_dp, 8.0_dp, 0.0_dp, 3.0_dp, 6.0_dp, &
      8.0_dp, 9.0_dp, 7.0_dp, 3.0_dp, 0.0_dp, 5.0_dp, &
      7.0_dp, 8.0_dp, 9.0_dp, 6.0_dp, 5.0_dp, 0.0_dp], [6, 6])

   call fastme_bal(distance, tree, info)
   if (info /= 0) error stop 'FastME reconstruction failed'

   print '(a)', 'parent child length'
   do e = 1, tree%nedge()
      print '(i0,1x,i0,1x,f10.5)', tree%edge(e, 1), tree%edge(e, 2), tree%edge_length(e)
   end do
end program fastme_example
