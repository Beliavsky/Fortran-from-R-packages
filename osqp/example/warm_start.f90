program warm_start
   use osqp
   implicit none
   real(dp) :: p(2,2), q(2), a(2,2), l(2), u(2)
   type(osqp_model) :: model
   type(osqp_solver) :: solver
   type(osqp_solution) :: cold, warm
   type(osqp_settings) :: settings
   integer(osqp_int) :: status

   if (.not. osqp_backend_available()) then
      print '(a)', 'Build the backend first with scripts/build_backend.bat or .sh.'
      stop
   end if
   p = reshape([2.0_dp,0.0_dp,0.0_dp,2.0_dp],[2,2])
   q = [-2.0_dp,-5.0_dp]
   a = reshape([1.0_dp,0.0_dp,0.0_dp,1.0_dp],[2,2])
   l = [0.0_dp,0.0_dp]
   u = [huge(1.0_dp),huge(1.0_dp)]
   settings%verbose = .false.
   call osqp_model_from_dense(model,q,status,p,a,l,u)
   call osqp_setup(solver,model,status,settings)
   call osqp_solve_solver(solver,cold,status)
   call osqp_warm_start(solver,status,x=cold%x,y=cold%y)
   call osqp_solve_solver(solver,warm,status)
   print '(a,i0)', 'cold iterations: ', cold%info%iter
   print '(a,i0)', 'warm iterations: ', warm%info%iter
   call osqp_cleanup(solver)
end program warm_start
