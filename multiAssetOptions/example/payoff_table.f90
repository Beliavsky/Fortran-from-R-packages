program payoff_table
   use multi_asset_options
   implicit none

   type(grid_set) :: grid
   type(status_type) :: status
   real(dp), allocatable :: values(:)
   integer :: i

   allocate(grid%asset(2),grid%dims(2),grid%strides(2))
   grid%asset(1)%x = [80.0_dp,100.0_dp,120.0_dp]
   grid%asset(2)%x = [80.0_dp,100.0_dp,120.0_dp]
   grid%dims = [3,3]
   grid%strides = [1,3]
   grid%n_nodes = 9

   call payoff_values(payoff_best_of,[option_call,option_call], &
      [100.0_dp,100.0_dp],grid,values,status)
   if (.not. status%ok()) error stop status%message

   print '(a)', 'Classically ordered best-of call payoff:'
   do i = 1, size(values)
      print '(i4,1x,f10.4)', i, values(i)
   end do
end program payoff_table
