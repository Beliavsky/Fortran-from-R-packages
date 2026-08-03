program persistent_solver_example
   use highs
   implicit none
   type(highs_model) :: model
   type(highs_solver) :: solver
   type(highs_solution) :: solution
   integer(highs_int) :: status
   real(dp) :: a(1,2)

   if (.not. highs_backend_available()) then
      print '(a)', "HiGHS backend not found. Run scripts/build_backend first."
      stop
   end if

   a = reshape([1.0_dp,1.0_dp], [1,2])
   call highs_model_from_dense(model, [1.0_dp,1.0_dp], [0.0_dp,0.0_dp], &
      [highs_default_infinity,highs_default_infinity], status, a=a, lhs=[1.0_dp], &
      rhs=[highs_default_infinity])
   call highs_new_solver(solver, status)
   call highs_pass_model(solver, model, status)
   call highs_run(solver, status)
   call highs_get_solution(solver, solution)
   print '(a,*(f8.4,1x))', "initial x: ", solution%col_value

   call highs_change_costs(solver, [1,2], [2.0_dp,1.0_dp], status)
   call highs_set_start(solver, solution%col_value, status)
   call highs_run(solver, status)
   call highs_get_solution(solver, solution)
   print '(a,*(f8.4,1x))', "after cost update: ", solution%col_value
   call highs_destroy_solver(solver)
end program persistent_solver_example
