program chrono_mpl_example
   use ape, only : chrono_mpl, dp, make_phylo_tree, phylo_tree
   implicit none

   integer :: edge(6, 2)
   real(dp) :: edge_length(6)
   real(dp), allocatable :: p_value(:)
   real(dp), allocatable :: standard_error(:)
   type(phylo_tree) :: dated_tree
   type(phylo_tree) :: tree
   integer :: info
   integer :: e

   edge = reshape([5, 5, 6, 6, 7, 7, 6, 7, 1, 2, 3, 4], [6, 2])
   edge_length = [4.0_dp, 6.0_dp, 1.0_dp, 3.0_dp, 2.0_dp, 2.0_dp]
   tree = make_phylo_tree(4, edge, edge_length)

   call chrono_mpl(tree, dated_tree, info, standard_error, p_value)
   if (info /= 0) error stop 'chronoMPL dating failed'

   print '(a)', 'dated parent child length'
   do e = 1, dated_tree%nedge()
      print '(i0,1x,i0,1x,f10.5)', dated_tree%edge(e, 1), dated_tree%edge(e, 2), dated_tree%edge_length(e)
   end do
   print '(a,*(1x,f9.5))', 'internal-node standard errors:', standard_error
   print '(a,*(1x,f9.5))', 'clock-test p-values:', p_value
end program chrono_mpl_example
