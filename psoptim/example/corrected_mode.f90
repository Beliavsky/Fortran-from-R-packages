program corrected_mode_example
   use psoptim, only : dp, ps_control, ps_result, ps_optimize
   implicit none

   type(ps_control) :: control
   type(ps_result) :: result
   real(dp) :: xmin(3), xmax(3), vmax(3)

   xmin = [-1.0_dp, 10.0_dp, -100.0_dp]
   xmax = [ 1.0_dp, 12.0_dp,  -90.0_dp]
   vmax = [0.5_dp, 0.5_dp, 2.0_dp]

   control%n = 60
   control%max_loop = 150
   control%legacy_initial_bounds = .false.
   control%legacy_velocity_initialization = .false.
   control%legacy_velocity_clip = .false.

   call ps_optimize(objective, xmin, xmax, vmax, result, control)
   print '(a,3f12.5)', 'solution: ', result%sol
   print '(a,es14.6)', 'fitness:  ', result%val

contains

   function objective(x) result(f)
      real(dp), intent(in) :: x(:)
      real(dp) :: f
      f = -((x(1)-0.25_dp)**2 + (x(2)-11.0_dp)**2 + (x(3)+95.0_dp)**2)
   end function objective

end program corrected_mode_example
