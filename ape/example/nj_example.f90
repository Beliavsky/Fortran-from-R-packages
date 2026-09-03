program nj_example
   use ape
   implicit none
   real(dp) :: distance(5, 5)
   type(phylo_tree) :: tree
   integer :: i
   integer :: info

   distance = reshape([ &
      0.0_dp, 5.0_dp, 9.0_dp, 9.0_dp, 8.0_dp, &
      5.0_dp, 0.0_dp, 10.0_dp, 10.0_dp, 9.0_dp, &
      9.0_dp, 10.0_dp, 0.0_dp, 8.0_dp, 7.0_dp, &
      9.0_dp, 10.0_dp, 8.0_dp, 0.0_dp, 3.0_dp, &
      8.0_dp, 9.0_dp, 7.0_dp, 3.0_dp, 0.0_dp], [5, 5])

   call nj(distance, tree, info)
   if (info /= 0) error stop 'neighbor joining failed'

   print '(a)', ' parent child       length'
   do i = 1, tree%nedge()
      print '(2i7,f13.6)', tree%edge(i, 1), tree%edge(i, 2), tree%edge_length(i)
   end do
end program nj_example
