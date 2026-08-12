program test_bounds_and_history
   use psoptim, only : dp, ps_control, ps_result, ps_optimize
   implicit none
   type(ps_control) :: control
   type(ps_result) :: result
   real(dp) :: xmin(3), xmax(3), vmax(3)
   integer :: j

   xmin = [-2.0_dp, 10.0_dp, -50.0_dp]
   xmax = [ 3.0_dp, 14.0_dp, -40.0_dp]
   vmax = [0.8_dp, 0.5_dp, 1.5_dp]
   control%n = 50
   control%max_loop = 120
   control%seed = 91
   control%legacy_initial_bounds = .false.
   control%legacy_velocity_initialization = .false.
   control%legacy_velocity_clip = .false.

   call ps_optimize(objective, xmin, xmax, vmax, result, control)
   do j = 1, 3
      if (minval(result%final_population(j,:)) < xmin(j)) error stop "lower bound violated"
      if (maxval(result%final_population(j,:)) > xmax(j)) error stop "upper bound violated"
      if (maxval(abs(result%final_velocity(j,:))) > vmax(j) + 1.0e-14_dp) &
         error stop "velocity limit violated"
   end do
   if (size(result%best_history) /= control%max_loop) error stop "bad best history size"
   if (size(result%mean_history) /= control%max_loop) error stop "bad mean history size"
   if (any(result%best_history(2:) < result%best_history(:control%max_loop-1))) &
      error stop "global-best history must be nondecreasing"
   if (result%val < -1.0e-4_dp) error stop "corrected-mode quadratic did not converge"
   print *, "test_bounds_and_history: PASS"
contains
   function objective(x) result(f)
      real(dp), intent(in) :: x(:)
      real(dp) :: f
      f = -((x(1)-0.5_dp)**2 + (x(2)-12.0_dp)**2 + (x(3)+45.0_dp)**2)
   end function objective
end program test_bounds_and_history
