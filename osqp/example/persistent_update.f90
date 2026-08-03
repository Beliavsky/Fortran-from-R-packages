program persistent_update
   use osqp
   implicit none
   real(dp) :: p(2,2), q(2), a(3,2), l(3), u(3)
   type(osqp_model) :: model
   type(osqp_solver) :: solver
   type(osqp_solution) :: solution
   type(osqp_settings) :: settings
   integer(osqp_int) :: status

   if (.not. osqp_backend_available()) then
      print '(a)', 'Build the backend first with scripts/build_backend.bat or .sh.'
      stop
   end if
   p = reshape([4.0_dp,1.0_dp,1.0_dp,2.0_dp],[2,2])
   q = [1.0_dp,1.0_dp]
   a = reshape([1.0_dp,1.0_dp,0.0_dp,1.0_dp,0.0_dp,1.0_dp],[3,2])
   l = [1.0_dp,0.0_dp,0.0_dp]
   u = [1.0_dp,0.7_dp,0.7_dp]
   settings%verbose = .false.
   call osqp_model_from_dense(model,q,status,p,a,l,u)
   call osqp_setup(solver,model,status,settings)
   call osqp_solve_solver(solver,solution,status)
   print '(a,*(f10.5,1x))', 'first x: ', solution%x
   call osqp_update(solver,status,q=[2.0_dp,3.0_dp])
   call osqp_solve_solver(solver,solution,status)
   print '(a,*(f10.5,1x))', 'updated x: ', solution%x
   call osqp_cleanup(solver)
end program persistent_update
