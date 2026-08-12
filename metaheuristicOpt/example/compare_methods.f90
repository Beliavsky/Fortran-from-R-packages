program compare_methods
   use metaheuristic_opt, only : dp, mh_control, mh_result, metaopt, sphere
   implicit none
   character(len=4), parameter :: names(4) = ['PSO ','GWO ','DE  ','WOA ']
   type(mh_control) :: control
   type(mh_result) :: result
   real(dp) :: lower(6), upper(6)
   integer :: i

   lower = -10.0_dp
   upper = 10.0_dp
   control%num_population = 40
   control%max_iter = 150
   control%legacy_quirks = .false.
   do i = 1, size(names)
      control%seed = 1000+i
      call metaopt(trim(names(i)), sphere, lower, upper, result, control)
      print '(a,2x,es14.6)', trim(names(i)), result%value
   end do
end program compare_methods
