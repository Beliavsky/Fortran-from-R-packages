program demo_highs
   use highs
   implicit none
   type(highs_solution) :: solution
   integer(highs_int) :: status
   real(dp) :: a(3,2), q(2,2)

   print '(a)', "highs-fortran demonstration"
   if (.not. highs_backend_available()) then
      print '(a)', "The Fortran frontend is installed, but the HiGHS backend library is absent."
      print '(a)', "Windows: scripts\\build_backend.bat"
      print '(a)', "Unix:   scripts/build_backend.sh"
      print '(a)', "Then run: fpm run"
      stop
   end if
   print '(a,a)', "HiGHS backend version: ", highs_backend_version()

   a = reshape([0.0_dp,1.0_dp,3.0_dp, 1.0_dp,2.0_dp,2.0_dp], [3,2])
   call highs_solve([1.0_dp,1.0_dp], [0.0_dp,1.0_dp], &
      [4.0_dp,highs_default_infinity], solution, status, a=a, &
      lhs=[-highs_default_infinity,5.0_dp,6.0_dp], &
      rhs=[7.0_dp,15.0_dp,highs_default_infinity], offset=3.0_dp)
   print '(a,*(f9.4,1x))', "LP x = ", solution%col_value
   print '(a,f10.4)', "LP objective = ", solution%objective_value

   q = reshape([1.0_dp,0.0_dp,0.0_dp,1.0_dp], [2,2])
   call highs_solve([-1.0_dp,-2.0_dp], [0.0_dp,0.0_dp], &
      [highs_default_infinity,highs_default_infinity], solution, status, q=q)
   print '(a,*(f9.4,1x))', "QP x = ", solution%col_value
   print '(a,f10.4)', "QP objective = ", solution%objective_value

   call highs_solve([3.0_dp,2.0_dp], [0.0_dp,0.0_dp], [1.0_dp,1.0_dp], &
      solution, status, a=reshape([2.0_dp,1.0_dp],[1,2]), &
      lhs=[-highs_default_infinity], rhs=[2.0_dp], &
      vartype=[highs_var_integer,highs_var_integer], maximum=.true.)
   print '(a,*(f6.1,1x))', "MIP x = ", solution%col_value
   print '(a,f8.2)', "MIP objective = ", solution%objective_value
end program demo_highs
