program node_spacing
   use multi_asset_options
   implicit none

   type(status_type) :: status
   real(dp), allocatable :: nodes(:)
   integer :: i

   call node_spacer(100.0_dp,0.0_dp,500.0_dp,26,5.0_dp,1,nodes,status)
   if (.not. status%ok()) error stop status%message

   print '(a)', 'Nonuniform finite-difference nodes:'
   do i = 1, size(nodes)
      print '(i4,1x,f12.5)', i, nodes(i)
   end do
end program node_spacing
