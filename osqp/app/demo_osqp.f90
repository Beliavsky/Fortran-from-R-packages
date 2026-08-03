program demo_osqp
   use osqp
   implicit none
   real(dp) :: p(2,2), q(2), a(3,2), l(3), u(3)
   type(osqp_solution) :: solution
   type(osqp_settings) :: settings
   integer(osqp_int) :: status

   if (.not. osqp_backend_available()) then
      print '(a)', 'The Fortran frontend is built, but the OSQP backend is not installed.'
      print '(a)', 'Run scripts\build_backend.bat on Windows or scripts/build_backend.sh on Unix.'
      print '(a)', 'Then run fpm run again.'
      stop
   end if
   print '(a,a)', 'OSQP backend version: ', osqp_backend_version()
   p = reshape([4.0_dp,1.0_dp,1.0_dp,2.0_dp],[2,2])
   q = [1.0_dp,1.0_dp]
   a = reshape([1.0_dp,1.0_dp,0.0_dp,1.0_dp,0.0_dp,1.0_dp],[3,2])
   l = [1.0_dp,0.0_dp,0.0_dp]
   u = [1.0_dp,0.7_dp,0.7_dp]
   settings%verbose = .false.
   settings%polishing = .true.
   call solve_osqp(q,solution,status,p=p,a=a,l=l,u=u,settings=settings)
   print '(a,a)', 'status: ', trim(solution%status)
   print '(a,*(f12.7,1x))', 'solution: ', solution%x
   print '(a,f12.7)', 'objective: ', solution%info%obj_val
end program demo_osqp
