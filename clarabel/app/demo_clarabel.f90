program demo_clarabel
   use clarabel
   implicit none
   real(dp) :: pd(2,2), ad(1,2), q(2), b(1)
   type(csc_matrix) :: p, a
   type(clarabel_cone) :: cones(1)
   type(clarabel_settings) :: settings
   type(clarabel_solution) :: sol
   integer :: code, env_status, env_length
   character(len=:), allocatable :: message

   ! Minimize 0.5*x'P*x + q'x subject to x1+x2=1.
   pd = reshape([4.0_dp, 1.0_dp, 1.0_dp, 2.0_dp], shape(pd))
   ad = reshape([1.0_dp, 1.0_dp], shape(ad))
   q = [1.0_dp, 1.0_dp]
   b = 1.0_dp
   p = csc_from_symmetric_upper(pd)
   a = csc_from_dense(ad)
   cones(1) = zero_cone(1)
   settings = default_clarabel_settings()
   settings%verbose = .false.
   call clarabel_solve_problem(p, q, a, b, cones, sol, settings, code, message)
   if (code /= 0) then
      call get_environment_variable("CLARABEL_FORTRAN_BRIDGE", length=env_length, status=env_status)
      if (code == clarabel_backend_unavailable .and. (env_status /= 0 .or. env_length == 0)) then
         print '(a)', 'Clarabel Fortran frontend: build successful.'
         print '(a)', 'The production Clarabel.rs backend DLL is not installed yet.'
         print '(a)', 'On Windows, run:'
         print '(a)', '  scripts\build_with_backend.bat run'
         print '(a)', 'or build it once and then use ordinary FPM commands:'
         print '(a)', '  scripts\build_backend.bat'
         print '(a)', '  fpm run'
         stop
      end if
      error stop message
   end if
   print '(a,a)', 'status: ', status_name(sol%status)
   print '(a,2f12.6)', 'x: ', sol%x
   print '(a,f12.6)', 'objective: ', sol%obj_val
   print '(a,es12.4)', 'primal residual: ', sol%r_prim
end program demo_clarabel
