program qp_example
   use highs
   implicit none
   real(dp) :: q(2,2)
   type(highs_solution) :: solution
   integer(highs_int) :: status

   if (.not. highs_backend_available()) then
      print '(a)', "HiGHS backend not found. Run scripts/build_backend first."
      stop
   end if
   q = reshape([1.0_dp,0.0_dp,0.0_dp,1.0_dp], [2,2])
   call highs_solve([-1.0_dp,-2.0_dp], [0.0_dp,0.0_dp], &
      [highs_default_infinity,highs_default_infinity], solution, status, q=q)
   print '(a,*(f10.5,1x))', "QP solution: ", solution%col_value
   print '(a,f12.6)', "objective: ", solution%objective_value
end program qp_example
