program test_legacy_semantics
   use psoptim, only : dp, ps_control, ps_result, ps_optimize
   implicit none
   type(ps_control) :: c_legacy, c_fixed
   type(ps_result) :: r_legacy, r_fixed
   real(dp) :: xmin(2), xmax(2), vmax(2)

   xmin = [0.0_dp, 100.0_dp]
   xmax = [1.0_dp, 101.0_dp]
   vmax = [0.01_dp, 0.01_dp]

   c_legacy%n = 20
   c_legacy%max_loop = 0
   c_legacy%seed = 12
   call ps_optimize(objective, xmin, xmax, vmax, r_legacy, c_legacy)
   if (maxval(r_legacy%final_population(2,:)) > 1.0_dp) &
      error stop "legacy initialization did not use dimension-1 bounds"

   c_fixed = c_legacy
   c_fixed%legacy_initial_bounds = .false.
   c_fixed%legacy_velocity_initialization = .false.
   c_fixed%legacy_velocity_clip = .false.
   call ps_optimize(objective, xmin, xmax, vmax, r_fixed, c_fixed)
   if (minval(r_fixed%final_population(2,:)) < 100.0_dp) &
      error stop "corrected initialization ignored dimension-2 bounds"
   if (maxval(r_fixed%final_population(2,:)) > 101.0_dp) &
      error stop "corrected initialization exceeds dimension-2 bounds"

   print *, "test_legacy_semantics: PASS"
contains
   function objective(x) result(f)
      real(dp), intent(in) :: x(:)
      real(dp) :: f
      f = -sum(x*x)
   end function objective
end program test_legacy_semantics
