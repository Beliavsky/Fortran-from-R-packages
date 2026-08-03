program test_basic_qp
   use osqp
   implicit none
   real(dp) :: p(2,2), q(2), a(5,2), l(5), u(5)
   type(osqp_solution) :: sol
   type(osqp_settings) :: settings
   integer(osqp_int) :: status

   if (.not. osqp_backend_available()) then
      print *, "SKIP test_basic_qp: ", trim(osqp_backend_error())
      stop
   end if
   p = reshape([11.0_dp,0.0_dp,0.0_dp,0.0_dp],[2,2])
   q = [3.0_dp,4.0_dp]
   a = reshape([-1.0_dp,0.0_dp,-1.0_dp,2.0_dp,3.0_dp, &
                0.0_dp,-1.0_dp,-3.0_dp,5.0_dp,4.0_dp],[5,2])
   l = -huge(1.0_dp)
   u = [0.0_dp,0.0_dp,-15.0_dp,100.0_dp,80.0_dp]
   settings%verbose = .false.
   settings%eps_abs = 1.0e-6_dp
   settings%eps_rel = 1.0e-6_dp
   call solve_osqp(q,sol,status,p=p,a=a,l=l,u=u,settings=settings)
   call check(status == 0, "solver call")
   call check(sol%solved(), "solved status")
   call check(maxval(abs(sol%x-[0.0_dp,5.0_dp])) < 2.0e-4_dp, "primal solution")
   call check(maxval(abs(sol%y-[1.6666667_dp,0.0_dp,1.3333333_dp,0.0_dp,0.0_dp])) < 3.0e-3_dp, "dual solution")
   call check(abs(sol%info%obj_val-20.0_dp) < 2.0e-3_dp, "objective")
   print *, "PASS test_basic_qp"
contains
   subroutine check(ok, message)
      logical, intent(in) :: ok
      character(len=*), intent(in) :: message
      if (.not. ok) error stop message
   end subroutine check
end program test_basic_qp
