program infeasibility
   use osqp
   implicit none
   real(dp) :: p(1,1), q(1), a(2,1), l(2), u(2)
   type(osqp_solution) :: solution
   type(osqp_settings) :: settings
   integer(osqp_int) :: status

   if (.not. osqp_backend_available()) then
      print '(a)', 'Build the backend first with scripts/build_backend.bat or .sh.'
      stop
   end if
   p=1.0_dp; q=0.0_dp; a=1.0_dp
   l=[1.0_dp,-huge(1.0_dp)]; u=[huge(1.0_dp),0.0_dp]
   settings%verbose=.false.
   call solve_osqp(q,solution,status,p=p,a=a,l=l,u=u,settings=settings)
   print '(a,a)', 'status: ', trim(solution%status)
   print '(a,*(es12.4,1x))', 'primal infeasibility certificate: ', solution%prim_inf_cert
end program infeasibility
