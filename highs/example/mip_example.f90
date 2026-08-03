program mip_example
   use highs
   implicit none
   real(dp) :: a(2,4)
   type(highs_solution) :: solution
   integer(highs_int) :: status

   if (.not. highs_backend_available()) then
      print '(a)', "HiGHS backend not found. Run scripts/build_backend first."
      stop
   end if

   ! Binary knapsack: maximize 8x1 + 6x2 + 5x3 + 4x4.
   a = reshape([4.0_dp,1.0_dp, 3.0_dp,1.0_dp, 2.0_dp,0.0_dp, 2.0_dp,1.0_dp], [2,4])
   call highs_solve([8.0_dp,6.0_dp,5.0_dp,4.0_dp], [0.0_dp,0.0_dp,0.0_dp,0.0_dp], &
      [1.0_dp,1.0_dp,1.0_dp,1.0_dp], solution, status, a=a, &
      lhs=[-highs_default_infinity,-highs_default_infinity], rhs=[6.0_dp,2.0_dp], &
      vartype=[highs_var_integer,highs_var_integer,highs_var_integer,highs_var_integer], &
      maximum=.true.)
   print '(a,*(f6.1,1x))', "selected: ", solution%col_value
   print '(a,f10.3)', "value: ", solution%objective_value
end program mip_example
