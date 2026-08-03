program test_infeasible
   use osqp
   implicit none
   real(dp) :: p(1,1), q(1), a(2,1), l(2), u(2)
   type(osqp_solution) :: sol
   type(osqp_settings) :: settings
   integer(osqp_int) :: status

   if (.not. osqp_backend_available()) then
      print *, "SKIP test_infeasible: ", trim(osqp_backend_error())
      stop
   end if
   p(1,1)=1.0_dp; q=0.0_dp
   a(:,1)=[1.0_dp,1.0_dp]
   l=[1.0_dp,-huge(1.0_dp)]
   u=[huge(1.0_dp),0.0_dp]
   settings%verbose=.false.; settings%eps_abs=1.0e-6_dp; settings%eps_rel=1.0e-6_dp
   call solve_osqp(q,sol,status,p=p,a=a,l=l,u=u,settings=settings)
   call check(status==0,"solve call")
   call check(sol%info%status_val==osqp_primal_infeasible .or. &
      sol%info%status_val==osqp_primal_infeasible_inaccurate,"primal infeasible status")
   call check(size(sol%prim_inf_cert)==2,"certificate dimension")
   print *, "PASS test_infeasible"
contains
   subroutine check(ok, message)
      logical, intent(in) :: ok
      character(len=*), intent(in) :: message
      if (.not. ok) error stop message
   end subroutine check
end program test_infeasible
