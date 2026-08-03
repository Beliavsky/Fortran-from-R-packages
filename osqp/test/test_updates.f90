program test_updates
   use osqp
   implicit none
   real(dp) :: p(2,2), q(2), a(5,2), l(5), u(5)
   real(dp) :: qnew(2), lnew(5), unew(5), px(2)
   type(osqp_model) :: model
   type(osqp_solver) :: solver
   type(osqp_solution) :: sol
   type(osqp_settings) :: settings
   integer(osqp_int) :: status

   if (.not. osqp_backend_available()) then
      print *, "SKIP test_updates: ", trim(osqp_backend_error())
      stop
   end if
   p = reshape([11.0_dp,0.0_dp,0.0_dp,1.0e-5_dp],[2,2])
   q = [3.0_dp,4.0_dp]
   a = reshape([-1.0_dp,0.0_dp,-1.0_dp,2.0_dp,3.0_dp, &
                0.0_dp,-1.0_dp,-3.0_dp,5.0_dp,4.0_dp],[5,2])
   l = -huge(1.0_dp)
   u = [0.0_dp,0.0_dp,-15.0_dp,100.0_dp,80.0_dp]
   settings%verbose=.false.; settings%eps_abs=1.0e-8_dp; settings%eps_rel=1.0e-8_dp
   call osqp_model_from_dense(model,q,status,p,a,l,u)
   call osqp_setup(solver,model,status,settings)
   call check(status==0,"setup")
   call osqp_solve_solver(solver,sol,status)
   call check(maxval(abs(sol%x-[0.0_dp,5.0_dp]))<1.0e-3_dp,"initial solve")

   qnew=[10.0_dp,20.0_dp]
   call osqp_update(solver,status,q=qnew)
   call check(status==0,"q update")
   call osqp_solve_solver(solver,sol,status)
   call check(abs(sol%info%obj_val-100.0_dp)<1.0e-2_dp,"updated q objective")

   call osqp_update(solver,status,q=q)
   call check(status==0,"restore q")
   lnew=-100.0_dp; unew=1000.0_dp
   call osqp_update(solver,status,l=lnew,u=unew)
   call check(status==0,"bounds update")
   call osqp_solve_solver(solver,sol,status)
   call check(maxval(abs(sol%x-[-0.12727273_dp,-19.94909091_dp]))<3.0e-3_dp,"updated bounds solution")

   call osqp_update(solver,status,q=q, l=l, u=u)
   px=[1.0_dp,3.0_dp]
   call osqp_update(solver,status,px=px)
   call check(status==0,"P update")
   call osqp_solve_solver(solver,sol,status)
   call check(maxval(abs(sol%x-[2.5_dp,4.1666667_dp]))<3.0e-3_dp,"updated P solution")
   call osqp_cleanup(solver)
   print *, "PASS test_updates"
contains
   subroutine check(ok, message)
      logical, intent(in) :: ok
      character(len=*), intent(in) :: message
      if (.not. ok) error stop message
   end subroutine check
end program test_updates
