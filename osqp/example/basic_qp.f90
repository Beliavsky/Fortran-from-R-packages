program basic_qp
   use osqp
   implicit none
   real(dp) :: p(2,2), q(2), a(5,2), l(5), u(5)
   type(osqp_solution) :: solution
   type(osqp_settings) :: settings
   integer(osqp_int) :: status

   if (.not. osqp_backend_available()) then
      print '(a)', 'OSQP backend unavailable. Run scripts/build_backend.bat or scripts/build_backend.sh.'
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
   call solve_osqp(q, solution, status, p=p, a=a, l=l, u=u, settings=settings)
   print '(a,a)', 'status: ', trim(solution%status)
   print '(a,*(f12.6,1x))', 'x: ', solution%x
   print '(a,f12.6)', 'objective: ', solution%info%obj_val
end program basic_qp
