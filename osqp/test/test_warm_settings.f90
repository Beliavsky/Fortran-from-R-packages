program test_warm_settings
   use osqp
   implicit none
   real(dp) :: p(2,2), q(2), a(5,2), l(5), u(5)
   type(osqp_model) :: model
   type(osqp_solver) :: solver
   type(osqp_solution) :: s1, s2
   type(osqp_settings) :: settings, got
   integer(osqp_int) :: status, n, m, cold_iter

   if (.not. osqp_backend_available()) then
      print *, "SKIP test_warm_settings: ", trim(osqp_backend_error())
      stop
   end if
   p=reshape([11.0_dp,0.0_dp,0.0_dp,0.0_dp],[2,2]); q=[3.0_dp,4.0_dp]
   a=reshape([-1.0_dp,0.0_dp,-1.0_dp,2.0_dp,3.0_dp,0.0_dp,-1.0_dp,-3.0_dp,5.0_dp,4.0_dp],[5,2])
   l=-huge(1.0_dp); u=[0.0_dp,0.0_dp,-15.0_dp,100.0_dp,80.0_dp]
   settings%verbose=.false.; settings%eps_abs=1.0e-6_dp; settings%eps_rel=1.0e-6_dp
   call osqp_model_from_dense(model,q,status,p,a,l,u)
   call osqp_setup(solver,model,status,settings)
   call osqp_solve_solver(solver,s1,status)
   cold_iter=s1%info%iter
   call osqp_warm_start(solver,status,x=s1%x,y=s1%y)
   call osqp_solve_solver(solver,s2,status)
   call check(s2%info%iter<=cold_iter,"warm start iterations")
   call check(maxval(abs(s2%x-s1%x))<1.0e-5_dp,"warm start solution")
   call osqp_cold_start(solver,status)
   call check(status==0,"cold start")
   settings%max_iter=2000
   settings%rho=0.2_dp
   call osqp_update_settings(solver,settings,status)
   call check(status==0,"settings update")
   call osqp_get_settings(solver,got,status)
   call check(got%max_iter==2000 .and. abs(got%rho-0.2_dp)<1.0e-14_dp,"settings round trip")
   call osqp_get_dimensions(solver,n,m,status)
   call check(n==2 .and. m==5,"dimensions")
   call check(len_trim(osqp_backend_version())>0,"version")
   call check(iand(osqp_capabilities(),osqp_capability_direct_solver)/=0,"direct capability")
   call osqp_cleanup(solver)
   print *, "PASS test_warm_settings"
contains
   subroutine check(ok, message)
      logical, intent(in) :: ok
      character(len=*), intent(in) :: message
      if (.not. ok) error stop message
   end subroutine check
end program test_warm_settings
